"""
Individual step implementations for scRNA pipeline
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ..base import PipelineStep
from ...utils.command import CommandRunner, RScriptRunner

class SCRNAFeatureExtractionStep(PipelineStep):
    """Feature extraction step - calls external bash script"""
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        """
        Initialize feature extraction step
        
        Args:
            workdir: Working directory
            script_dir: Directory containing external scripts
            config: Configuration dictionary
        """
        super().__init__("01_features", workdir, config)
        self.script_dir = Path(script_dir) / "scrna" / "feature_extraction"
        self.main_script = self.script_dir / "scRNA_sites_grep_all_features_and_filtration_withASE_and_patched.whole_steps.sh"
        self.runner = CommandRunner(workdir=self.workdir)
        self.r_runner = RScriptRunner(workdir=self.workdir)
        
    def _execute(self, sample_id: str, mutation_list: Path, bam_file: Path,
                 barcode_file: Path, read_len: int = 91, 
                 running_type: str = "benchmark",
                 ase_filepath: Optional[str] = None,
                 threads: int = 4) -> Dict[str, Any]:
        """
        Execute feature extraction
        
        Returns:
            Dictionary with step results
        """
        self.logger.info(f"Extracting features for {sample_id}")
        
        # Check if script exists
        if not self.main_script.exists():
            raise FileNotFoundError(f"Script not found: {self.main_script}")
        
        # Prepare arguments matching the original script's 10 parameters
        args = [
            sample_id,
            str(mutation_list.absolute()),
            str(bam_file.absolute()),
            str(read_len),
            str(self.workdir.absolute()),           # out_dir_name
            str(bam_file.absolute()),               # filtered_bam
            str(threads),
            str(barcode_file.absolute()),
            ase_filepath if ase_filepath else "no",
            running_type
        ]
        
        self.logger.info(f"Running feature extraction script with {len(args)} parameters")
        self.runner.run_script(self.main_script, args)
        
        # Find output files matching real output structure
        patched_file = self.workdir / f"{sample_id}.{running_type}_patched.feature.txt"
        if patched_file.exists():
            self.outputs['features_file'] = patched_file
            self.logger.info(f"Found features file: {patched_file}")
        
        # Record depth_in_spots directory
        depth_dir = self.workdir / "depth_in_spots"
        if depth_dir.exists():
            self.outputs['depth_dir'] = depth_dir
        
        return {
            'sample_id': sample_id,
            'features_file': str(patched_file) if patched_file.exists() else None,
        }


class SCRNATreeInputStep(PipelineStep):
    """Tree input generation step"""
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        """
        Initialize tree input step
        
        Args:
            workdir: Working directory
            script_dir: Directory containing external scripts
            config: Configuration dictionary
        """
        super().__init__("02_treeinput", workdir, config)
        self.script_dir = Path(script_dir) / "scrna" / "tree_input"
        self.main_script = self.script_dir / "run_PhyloSOLID.treeinput_data.scRNAmode_using_identifiers.sh"
        self.runner = CommandRunner(workdir=self.workdir)
        
    def _execute(self, sample_id: str, mutation_list: Path, bam_file: Path,
                 barcode_file: Path, cellnum: int, threads: int = 2) -> Dict[str, Any]:
        """
        Execute tree input generation
        
        Returns:
            Dictionary with step results
        """
        self.logger.info(f"Generating tree input for {sample_id}")
        
        # Check if script exists
        if not self.main_script.exists():
            raise FileNotFoundError(f"Script not found: {self.main_script}")
        
        # Workpath should be parent directory
        workpath = self.workdir.parent.absolute()
        
        # Prepare arguments
        args = [
            str(workpath),
            sample_id,
            str(barcode_file.absolute()),
            str(bam_file.absolute()),
            str(mutation_list.absolute()),
            str(cellnum),
            str(threads)
        ]
        
        # Run script
        self.runner.run_script(self.main_script, args)
        
        # Collect output files matching real structure
        treeinput_dir = self.workdir / 'treeinput'
        outputs = {}
        
        # Spot matrix file
        spot_files = list(treeinput_dir.glob(f"treeinput_spot_c_{cellnum}.csv"))
        if spot_files:
            self.outputs['spot_matrix'] = spot_files[0]
            outputs['spot_matrix'] = str(spot_files[0])
        
        # SCID barcode file
        scid_file = treeinput_dir / "treeinput_scid_barcode.txt"
        if scid_file.exists():
            self.outputs['scid_barcode'] = scid_file
            outputs['scid_barcode'] = str(scid_file)
        
        # Features file
        features_file = treeinput_dir / "features_file.txt"
        if features_file.exists():
            self.outputs['features_file'] = features_file
            outputs['features_file'] = str(features_file)
        
        return {
            'sample_id': sample_id,
            'outputs': outputs
        }


class SCRNATreeBuildingStep(PipelineStep):
    """Tree building step"""
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        """
        Initialize tree building step
        
        Args:
            workdir: Working directory
            script_dir: Directory containing external scripts
            config: Configuration dictionary
        """
        super().__init__("03_tree_building", workdir, config)
        self.script_dir = Path(script_dir) / "scrna" / "tree_building"
        self.main_script = self.script_dir / "bash_run_other_methods_for_scRNA_part1.sh"
        self.runner = CommandRunner(workdir=self.workdir)
        self.r_runner = RScriptRunner(workdir=self.workdir)
        
    def _execute(self, sample_id: str, cellnum: int, 
                 celltype_file: Optional[Path] = None,
                 metadata_file: Optional[Path] = None,
                 barcode_file: Optional[Path] = None) -> Dict[str, Any]:
        """
        Execute tree building
        
        Returns:
            Dictionary with step results
        """
        self.logger.info(f"Building tree for {sample_id}")
        
        # Prepare cell type file if needed
        if celltype_file is None and metadata_file and barcode_file:
            celltype_file = self._prepare_celltype_file(sample_id, metadata_file, barcode_file)
        
        # Check if script exists
        if not self.main_script.exists():
            raise FileNotFoundError(f"Script not found: {self.main_script}")
        
        # Workpath should be parent directory
        workpath = self.workdir.parent.absolute()
        
        # Prepare arguments
        args = [
            str(workpath),
            sample_id,
            str(cellnum),
            str(celltype_file) if celltype_file else "None"
        ]
        
        # Run script
        self.runner.run_script(self.main_script, args)
        
        # Find output files matching real structure
        results_dir = self.workdir / 'results'
        tree_file = results_dir / 'PhyloSOLID' / 'celltree.newick'
        
        if tree_file.exists():
            self.outputs['tree_file'] = tree_file
            self.logger.info(f"Tree file generated: {tree_file}")
        
        cfmatrix = results_dir / 'PhyloSOLID' / 'cell_by_mut.CFMatrix'
        if cfmatrix.exists():
            self.outputs['cfmatrix'] = cfmatrix
        
        return {
            'sample_id': sample_id,
            'tree_file': str(tree_file) if tree_file.exists() else None,
            'results_dir': str(results_dir)
        }
    
    def _prepare_celltype_file(self, sample_id: str, metadata_file: Path, 
                               barcode_file: Path) -> Optional[Path]:
        """
        Prepare cell type file from metadata
        
        Args:
            sample_id: Sample identifier
            metadata_file: Path to metadata file
            barcode_file: Path to barcode file
            
        Returns:
            Path to generated cell type file or None
        """
        output_file = self.workdir / f"celltype_file_for_{sample_id}.txt"
        
        # Path to cell type extraction script
        script_path = Path(__file__).parent.parent.parent / 'scripts' / 'python' / 'catch_celltype_from_metadata.py'
        
        if script_path.exists():
            import subprocess
            import sys
            
            cmd = [
                sys.executable,
                str(script_path),
                '-i1', str(metadata_file),
                '-i2', str(barcode_file),
                '-p', sample_id,
                '-o', str(output_file)
            ]
            
            subprocess.run(cmd, check=True)
            return output_file if output_file.exists() else None
        
        self.logger.warning(f"Cell type extraction script not found: {script_path}")
        return None

### `pipelines/scdna/__init__.py`
```python
"""scDNA-seq pipeline module"""

from .pipeline import scDNAPipeline

__all__ = ['scDNAPipeline']