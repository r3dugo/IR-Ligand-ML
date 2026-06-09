#!/usr/bin/env python3
"""
Batch LoQI conformer generation for CSV files with columns: id, smiles.

This script calls a local LoQI checkout without modifying or vendoring LoQI.
It writes SDF conformers, JSON metadata, batch result summaries, logs, and
failed-molecule records. Optionally converts SDF outputs to XYZ with RDKit.

Assumes structure:
Generate_LoQI/
├── generate_loqi_generic.py
└── LoQI/
    ├── scripts/sample_conformers.py
    ├── scripts/conf/loqi/loqi.yaml
    └── data/loqi.ckpt
"""

import sys
import csv
from pathlib import Path
import json
import time
import subprocess
import tempfile
import logging
from datetime import datetime
import signal
import os
import argparse

# Set up logging
def setup_logging(output_dir):
    """Set up comprehensive logging"""
    log_file = output_dir / f"loqi_molecule_generation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return log_file

class LoQIBatchProcessor:
    def __init__(self, output_dir, n_conformers=12, batch_size=50, sample_script=None, config_path=None, ckpt_path=None, convert_sdf_to_xyz=False):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.n_conformers = n_conformers
        self.batch_size = batch_size

        # LoQI sample script/config/ckpt are passed from CLI/main()
        self.sample_script = Path(sample_script)
        self.config_path = Path(config_path)
        self.ckpt_path = Path(ckpt_path)

        self.convert_sdf_to_xyz = bool(convert_sdf_to_xyz)
        
        # Progress tracking
        self.processed = 0
        self.failed = 0
        self.total_time = 0
        self.start_time = None
        
        # Results storage
        self.results = []
        self.failed_molecules = []
        
        # Checkpointing
        self.checkpoint_file = self.output_dir / "checkpoint.json"
        self.results_file = self.output_dir / "batch_results.csv"
        self.failed_file = self.output_dir / "failed_molecules.csv"
        
        # Set up logging
        self.log_file = setup_logging(self.output_dir)
        self.logger = logging.getLogger(__name__)
        
        # Signal handling for graceful shutdown
        signal.signal(signal.SIGINT, self.signal_handler)
        signal.signal(signal.SIGTERM, self.signal_handler)
        
        self.logger.info(f"Initialised LoQI batch processor")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Conformers per molecule: {self.n_conformers}")
        self.logger.info(f"Convert SDF->XYZ: {self.convert_sdf_to_xyz}")
        self.logger.info(f"Log file: {self.log_file}")
    
    def signal_handler(self, signum, frame):
        """Handle interruption gracefully"""
        self.logger.warning(f"Received signal {signum}. Saving progress and exiting...")
        self.save_checkpoint()
        self.save_results()
        sys.exit(0)
    
    def load_checkpoint(self):
        """Load previous progress if available"""
        if self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file, 'r') as f:
                    checkpoint = json.load(f)
                
                self.processed = checkpoint.get('processed', 0)
                self.failed = checkpoint.get('failed', 0)
                self.total_time = checkpoint.get('total_time', 0)
                
                self.logger.info(f"Loaded checkpoint: {self.processed} processed, {self.failed} failed")
                return checkpoint.get('last_index', 0)
            except Exception as e:
                self.logger.warning(f"Failed to load checkpoint: {e}")
        
        return 0
    
    def save_checkpoint(self):
        """Save current progress"""
        checkpoint = {
            'processed': self.processed,
            'failed': self.failed,
            'total_time': self.total_time,
            'last_index': self.processed + self.failed,
            'timestamp': datetime.now().isoformat()
        }
        
        try:
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint, f, indent=2)
        except Exception as e:
            self.logger.error(f"Failed to save checkpoint: {e}")
    
    def save_results(self):
        """Save current results to files"""
        # Save successful results
        if self.results:
            try:
                fieldnames = ['id', 'smiles', 'n_conformers', 'energy_range_kcal', 'generation_time', 'sdf_file', 'xyz_file', 'json_file']
                with open(self.results_file, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    for r in self.results:
                        writer.writerow({k: r.get(k, '') for k in fieldnames})
                self.logger.info(f"Saved {len(self.results)} results to {self.results_file}")
            except Exception as e:
                self.logger.error(f"Failed to save results: {e}")
        
        # Save failed molecules
        if self.failed_molecules:
            try:
                fieldnames = ['id', 'smiles', 'error']
                with open(self.failed_file, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    for fm in self.failed_molecules:
                        if isinstance(fm, dict):
                            row = {k: fm.get(k, '') for k in fieldnames}
                        else:
                            # assume tuple (id, smiles, error)
                            row = {
                                'id': fm[0] if len(fm) > 0 else '',
                                'smiles': fm[1] if len(fm) > 1 else '',
                                'error': fm[2] if len(fm) > 2 else ''
                            }
                        writer.writerow(row)
                self.logger.info(f"Saved {len(self.failed_molecules)} failed molecules to {self.failed_file}")
            except Exception as e:
                self.logger.error(f"Failed to save failed molecules: {e}")
    
    def generate_conformers_direct(self, smiles, id):
        """Generate conformers by invoking LoQI's `sample_conformers.py` locally.

        The sample script writes SDF output; we parse the SDF and optionally convert
        to XYZ. Legacy remote client support has been removed.
        """
        suffix = '.sdf'

        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            tmp_path = tmp.name

        try:
            # Ensure sample script, config and checkpoint exist
            if not self.sample_script.exists():
                return {'success': False, 'error': f"LoQI sample script not found: {self.sample_script}", 'n_conformers': 0}
            if not self.config_path.exists():
                return {'success': False, 'error': f"LoQI config not found: {self.config_path}", 'n_conformers': 0}
            if not self.ckpt_path.exists():
                return {'success': False, 'error': f"LoQI checkpoint not found: {self.ckpt_path}", 'n_conformers': 0}

            cmd = [
                sys.executable,
                str(self.sample_script),
                '--config', str(self.config_path),
                '--ckpt', str(self.ckpt_path),
                '--input', smiles,
                '--output', tmp_path,
                '--n_confs', str(self.n_conformers),
                '--batch_size', '1'
            ]

            # Execute client in LoQI repo root so imports/relative paths resolve.
            try:
                loqi_root = self.sample_script.parents[1]
            except Exception:
                loqi_root = self.sample_script.parent

            env = os.environ.copy()
            src_path = str(loqi_root / 'src')
            existing_py = env.get('PYTHONPATH', '')
            env['PYTHONPATH'] = f"{src_path}:{existing_py}" if existing_py else src_path

            self.logger.debug(f"Running LoQI sample script in cwd={loqi_root} with PYTHONPATH={env['PYTHONPATH']}")

            # Execute client
            start_time = time.time()
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
                cwd=str(loqi_root),
                env=env,
            )
            generation_time = time.time() - start_time

            if result.returncode != 0:
                return {'success': False, 'error': (
                    f"Client error (code {result.returncode})\n"
                    f"STDERR:\n{result.stderr}\n"
                    f"STDOUT:\n{result.stdout}"
                ), 'n_conformers': 0, 'generation_time': generation_time}

            # Read output
            try:
                with open(tmp_path, 'r', encoding='utf-8', errors='replace') as f:
                    out_content = f.read()
            except Exception as e:
                return {'success': False, 'error': f"Failed to read output: {e}", 'n_conformers': 0, 'generation_time': generation_time}

            # Parse SDF output (sample_conformers.py writes SDF)
            parsed = self.parse_sdf_content(out_content, smiles, generation_time)

            if parsed.get('success'):
                self.save_individual_files(id, parsed)

            return parsed

        except subprocess.TimeoutExpired:
            return {'success': False, 'error': 'Timeout after 600 seconds', 'n_conformers': 0}
        except Exception as e:
            return {'success': False, 'error': f'Unexpected error: {str(e)}', 'n_conformers': 0}
        finally:
            try:
                Path(tmp_path).unlink()
            except Exception:
                pass
    
    # Legacy XYZ parsing removed: LoQI is invoked via sample_conformers.py which writes SDF.
    # If XYZ output is needed, use --convert-sdf-to-xyz to convert SDF -> XYZ via RDKit.

    def parse_sdf_content(self, sdf_content, smiles, generation_time):
        """Parse SDF content and extract conformers and energies.

        Looks for property blocks like:

        >  <Energy>
        -123.456

        Falls back to searching for any numeric value in property fields.
        """
        if not sdf_content.strip():
            return {'success': False, 'error': 'Empty SDF content', 'n_conformers': 0}

        # Split records by SDF delimiter
        records = [r.strip() for r in sdf_content.split('$$$$') if r.strip()]
        sdf_blocks = records
        energies = []

        import re
        prop_re = re.compile(r'>\s*<([^>]+)>\s*\n([^\n]+)', re.IGNORECASE)
        num_re = re.compile(r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?')

        for rec in records:
            energy_val = None
            # Only accept explicit property fields that mention 'energy' in their name.
            for m in prop_re.finditer(rec):
                prop_name = m.group(1).strip().lower()
                prop_val = m.group(2).strip()
                if 'energy' in prop_name:
                    num_m = num_re.search(prop_val)
                    if num_m:
                        try:
                            energy_val = float(num_m.group(0))
                        except Exception:
                            energy_val = None
                    break

            # Do NOT attempt to extract numbers from arbitrary parts of the SDF record
            # (coordinates, counts, bond indices, etc.). If no explicit energy property
            # is present, default to 0.0 to avoid inventing spurious energies.
            if energy_val is None:
                energy_val = 0.0

            energies.append(energy_val)

        # Convert to relative energies in kcal/mol if any non-zero found
        if energies and any(e != 0.0 for e in energies):
            min_e = min(energies)
            energies_kcal = [ (e - min_e) * 627.509 for e in energies ]
        else:
            energies_kcal = [0.0] * len(sdf_blocks)

        return {
            'success': True,
            'n_conformers': len(sdf_blocks),
            'smiles': smiles,
            'generation_time': generation_time,
            'sdf_data': sdf_content,
            'sdf_blocks': sdf_blocks,
            'energies': energies,
            'energies_kcal': energies_kcal
        }

    def convert_sdf_to_xyz_content(self, sdf_content, energies_kcal=None):
        """Convert SDF records to a multi-block XYZ string using RDKit.

        Returns a string containing concatenated XYZ blocks. Raises RuntimeError
        if RDKit is not available or parsing fails.
        """
        try:
            from rdkit import Chem
        except Exception as e:
            raise RuntimeError(
                "RDKit is required for SDF->XYZ conversion. Install via conda: 'conda install -c conda-forge rdkit'"
            ) from e

        records = [r.strip() for r in sdf_content.split('$$$$') if r.strip()]
        xyz_blocks = []

        for idx, rec in enumerate(records):
            try:
                mol = Chem.MolFromMolBlock(rec, sanitize=False, removeHs=False)
                if mol is None:
                    mol = Chem.MolFromMolBlock(rec, sanitize=True, removeHs=False)
                if mol is None:
                    raise ValueError('RDKit failed to parse SDF record')

                conf = mol.GetConformer()
                n_atoms = mol.GetNumAtoms()

                energy_str = ''
                if energies_kcal and idx < len(energies_kcal):
                    energy_str = f" energy_kcal={energies_kcal[idx]:.6f}"

                comment = f"Conformer{energy_str}"
                lines = [str(n_atoms), comment]

                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    symbol = atom.GetSymbol()
                    lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

                xyz_blocks.append('\n'.join(lines))
            except Exception as e:
                self.logger.warning(f"SDF->XYZ conversion: failed record {idx}: {e}")
                continue

        xyz_data = '\n'.join(block + '\n' for block in xyz_blocks)
        return xyz_data
    
    # extract_energy_from_comment removed: energies are read from SDF property fields only.
    
    def save_individual_files(self, id, result):
        """Save individual SDF/XYZ and JSON metadata files.

        If `self.convert_sdf_to_xyz` is True and an SDF is present, attempt
        an RDKit-based conversion to XYZ and save both files.
        """
        sdf_filename = ''
        xyz_filename = ''

        try:
            # Save SDF (if present)
            if 'sdf_data' in result:
                sdf_file = self.output_dir / f"{id}_conformers.sdf"
                sdf_file.write_text(result['sdf_data'])
                sdf_filename = sdf_file.name

                # Optionally convert to XYZ
                if self.convert_sdf_to_xyz:
                    try:
                        xyz_data = self.convert_sdf_to_xyz_content(result['sdf_data'], result.get('energies_kcal'))
                        xyz_file = self.output_dir / f"{id}_conformers.xyz"
                        xyz_file.write_text(xyz_data)
                        xyz_filename = xyz_file.name
                        result['xyz_data'] = xyz_data
                    except Exception as e:
                        self.logger.warning(f"Failed to convert SDF->XYZ for {id}: {e}")

            # If no SDF but we have XYZ data (legacy client), save XYZ
            elif 'xyz_data' in result:
                xyz_file = self.output_dir / f"{id}_conformers.xyz"
                xyz_file.write_text(result.get('xyz_data', ''))
                xyz_filename = xyz_file.name

            # Save JSON metadata (compact)
            json_file = self.output_dir / f"{id}_conformers.json"
            energies = result.get('energies_kcal', []) or []
            energy_range = (max(energies) - min(energies)) if energies else 0
            json_data = {
                'id': id,
                'smiles': result.get('smiles'),
                'n_conformers': result.get('n_conformers', 0),
                'energies_kcal': energies,
                'generation_time': result.get('generation_time', 0),
                'energy_range': energy_range,
                'sdf_file': sdf_filename,
                'xyz_file': xyz_filename
            }

            with open(json_file, 'w') as f:
                json.dump(json_data, f, indent=1)

            # Update result dict with filenames for downstream records
            result['sdf_file'] = sdf_filename
            result['xyz_file'] = xyz_filename
            result['json_file'] = json_file.name

            return result

        except Exception as e:
            self.logger.warning(f"Failed to save individual files for {id}: {e}")
            return result
    
    def process_dataset(self, csv_file):
        """Process the full dataset"""
        
        # Ensure the local LoQI sample script and its config/checkpoint are available
        if not self.sample_script.exists():
            self.logger.error(f"LoQI sample script not found: {self.sample_script}")
            return False
        self.logger.info(f"Using local LoQI sample script: {self.sample_script}")
        if not self.config_path.exists():
            self.logger.error(f"LoQI config not found: {self.config_path}")
            return False
        if not self.ckpt_path.exists():
            self.logger.error(f"LoQI checkpoint not found: {self.ckpt_path}")
            return False
        
        # Load dataset using csv to avoid pandas dependency
        self.logger.info(f"Loading molecule dataset from {csv_file}")
        rows = []
        try:
            with open(csv_file, newline='', encoding='utf-8', errors='replace') as f:
                reader = csv.DictReader(f)
                for r in reader:
                    # Normalise keys to lowercase for consistent access
                    r_norm = { (k.strip().lower() if k else k): v for k, v in r.items() }
                    rows.append(r_norm)
        except Exception as e:
            self.logger.error(f"Failed to read CSV {csv_file}: {e}")
            return False

        # Check required CSV columns
        if not rows:
            self.logger.error(f"CSV is empty: {csv_file}")
            return False

        required_cols = {"id", "smiles"}
        missing_cols = required_cols - set(rows[0].keys())

        if missing_cols:
            self.logger.error(f"CSV is missing required columns: {sorted(missing_cols)}")
            return False

        total_molecules = len(rows)
        self.logger.info(f"Loaded {total_molecules} molecules")

        # Load checkpoint
        start_index = self.load_checkpoint()
        if start_index > 0:
            self.logger.info(f"Resuming from index {start_index}")

        # Start processing
        self.start_time = time.time()

        for idx, row in enumerate(rows):
            # Skip already processed
            if idx < start_index:
                continue

            id = row.get('id')
            smiles = row.get('smiles')
            
            # Progress update
            progress = (idx + 1) / total_molecules * 100
            elapsed = time.time() - self.start_time
            eta = (elapsed / (idx + 1 - start_index)) * (total_molecules - idx - 1) if idx > start_index else 0
            
            self.logger.info(f"[{idx+1}/{total_molecules}] ({progress:.1f}%) Processing {id}: {smiles}")
            self.logger.info(f"   Elapsed: {elapsed/60:.1f}min, ETA: {eta/60:.1f}min")
            
            try:
                # Generate conformers
                result = self.generate_conformers_direct(smiles, id)
                
                if result['success']:
                    self.processed += 1
                    self.total_time += result['generation_time']
                    
                    energy_range = max(result['energies_kcal']) - min(result['energies_kcal']) if result['energies_kcal'] else 0
                    
                    self.logger.info(f"   Generated {result['n_conformers']} conformers ({result['generation_time']:.2f}s, {energy_range:.2f} kcal/mol range)")
                    
                    # Add to results
                    self.results.append({
                        'id': id,
                        'smiles': smiles,
                        'n_conformers': result['n_conformers'],
                        'energy_range_kcal': energy_range,
                        'generation_time': result['generation_time'],
                        'sdf_file': result.get('sdf_file', ''),
                        'xyz_file': result.get('xyz_file', ''),
                        'json_file': result.get('json_file', f"{id}_conformers.json")
                    })
                    
                else:
                    self.failed += 1
                    error_msg = result.get('error', 'Unknown error')
                    self.logger.warning(f"   Failed: {error_msg}")
                    self.failed_molecules.append((id, smiles, error_msg))
                
            except Exception as e:
                self.failed += 1
                error_msg = str(e)[:100]
                self.logger.error(f"    Exception: {error_msg}")
                self.failed_molecules.append((id, smiles, error_msg))
            
            # Save progress periodically
            if (idx + 1) % self.batch_size == 0:
                self.save_checkpoint()
                self.save_results()
                self.logger.info(f"   Checkpoint saved (processed: {self.processed}, failed: {self.failed})")
        
        # Final save
        self.save_checkpoint()
        self.save_results()
        
        # Summary
        total_time = time.time() - self.start_time
        self.print_summary(total_molecules, total_time)
        
        return self.processed > 0
    
    def print_summary(self, total_molecules, total_time):
        """Print final summary"""
        total_attempted = self.processed + self.failed
        success_rate = self.processed / total_attempted * 100 if total_attempted > 0 else 0.0

        self.logger.info("\n" + "="*80)
        self.logger.info("LoQI BATCH CONFORMER GENERATION COMPLETE")
        self.logger.info("="*80)
        self.logger.info(f"Successfully processed: {self.processed} molecules")
        self.logger.info(f"Failed: {self.failed} molecules")
        self.logger.info(f"Success rate: {success_rate:.1f}%")
        self.logger.info(f"Total wall time: {total_time/3600:.1f} hours")
        self.logger.info(f"Average generation time: {self.total_time/self.processed:.2f}s per molecule" if self.processed > 0 else "")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Results file: {self.results_file}")
        self.logger.info(f"Failed molecules: {self.failed_file}")
        
        if self.processed > 0:
            total_conformers = sum(r['n_conformers'] for r in self.results)
            avg_conformers = total_conformers / self.processed
            avg_energy_range = sum(r['energy_range_kcal'] for r in self.results) / self.processed
            
            self.logger.info(f"\nStatistics:")
            self.logger.info(f"   Total conformers generated: {total_conformers:,}")
            self.logger.info(f"   Average conformers per molecule: {avg_conformers:.1f}")
            self.logger.info(f"   Average energy range: {avg_energy_range:.2f} kcal/mol")
        
        self.logger.info("="*80)

def main():

    print("LoQI Large-Scale Molecule Conformer Generation")
    print("="*50)

    # CLI: allow choosing CSVs and output location
    parser = argparse.ArgumentParser(description="Run LoQI conformer generation on one or more CSVs")

    parser.add_argument('--csvs', nargs='+', required=True,
                        help='One or more CSV filenames or paths to process')
    
    parser.add_argument('--output-base', default=None,
                        help='Base output directory. Per-file dirs will be created inside this. Defaults to script dir')
    
    parser.add_argument('--n-conformers', type=int, default=12, help='Number of conformers to generate')

    parser.add_argument('--batch-size', type=int, default=50, help='Checkpoint batch size')

    # Scripts should be found in the LoQI repo relative to this script,
    # but allow overrides via CLI for flexibility.
    script_dir = Path(__file__).resolve().parent
    default_loqi_dir = script_dir / "LoQI"

    parser.add_argument(
        "--sample-script",
        default=str(default_loqi_dir / "scripts" / "sample_conformers.py"),
        help="Path to LoQI sample_conformers.py"
    )

    parser.add_argument(
        "--config",
        default=str(default_loqi_dir / "scripts" / "conf" / "loqi" / "loqi.yaml"),
        help="Path to LoQI config YAML"
    )

    parser.add_argument(
        "--ckpt",
        default=str(default_loqi_dir / "data" / "loqi.ckpt"),
        help="Path to LoQI checkpoint"
    )

    parser.add_argument(
        "--convert-sdf-to-xyz",
        action='store_true',
        help="If set, convert SDF outputs to multi-block XYZ files (requires RDKit)"
    )

    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent

    # Check LoQI files before starting any processing
    required_files = {
        "LoQI sample script": Path(args.sample_script),
        "LoQI config": Path(args.config),
        "LoQI checkpoint": Path(args.ckpt),
    }

    missing_files = [
        f"{label}: {path}"
        for label, path in required_files.items()
        if not path.exists()
    ]

    if missing_files:
        print("Missing required LoQI files:")
        for item in missing_files:
            print(f"  - {item}")
        print("\nEither clone/download LoQI into:")
        print(f"  {default_loqi_dir}")
        print("\nor pass paths manually with:")
        print("  --sample-script /path/to/sample_conformers.py")
        print("  --config /path/to/loqi.yaml")
        print("  --ckpt /path/to/loqi.ckpt")
        return

    # Resolve and check CSVs before running LoQIBatchProcessor
    csv_paths = []

    for csv_entry in args.csvs:
        csv_path = Path(csv_entry)
        if not csv_path.is_absolute():
            csv_path = base_dir / csv_path

        if not csv_path.exists():
            print(f"Input file not found: {csv_path}")
            continue

        csv_paths.append(csv_path)

    if not csv_paths:
        print("No valid input CSV files found. Exiting.")
        return

    # Resolve output base
    if args.output_base:
        out_base = Path(args.output_base)
    else:
        out_base = base_dir

    out_base.mkdir(parents=True, exist_ok=True)

    for csv_path in csv_paths:
        output_dir = Path(out_base) / f"conformers_full_loqi_{args.n_conformers}conf_{csv_path.stem}"

        # Initialise processor for this dataset
        processor = LoQIBatchProcessor(
            output_dir=output_dir,
            n_conformers=args.n_conformers,
            batch_size=args.batch_size,
            sample_script=args.sample_script,
            config_path=args.config,
            ckpt_path=args.ckpt,
            convert_sdf_to_xyz=args.convert_sdf_to_xyz,
        )

        print(f"\nProcessing {csv_path} -> output: {output_dir}")
        success = processor.process_dataset(str(csv_path))

        if success:
            print(f"\nCompleted: {csv_path.name} (output: {output_dir})")
        else:
            print(f"\nFailed: {csv_path.name} (see log: {processor.log_file})")

if __name__ == "__main__":
    main()
