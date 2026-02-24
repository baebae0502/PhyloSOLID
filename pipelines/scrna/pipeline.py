"""
Complete scRNA pipeline orchestrating all steps
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional, List
from concurrent.futures import ThreadPoolExecutor

from ..base import Pipeline
from .steps import (
    SCRNAFeatureExtractionStep,
    SCRNATreeInputStep,
    SCRNATreeBuildingStep
)

class SCRNAPipeline(Pipeline):
    """Complete scRNA pipeline"""
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        """
        Initialize scRNA pipeline
        
        Args:
            workdir: Working directory
            script_dir: Directory containing external scripts
            config: Configuration dictionary
        """
        super().__init__(workdir, config)
        self.script_dir = script_dir
        self.logger = logging.getLogger(__name__)
        
        # Initialize steps with numeric prefixes to indicate order
        self.add_step('feature_extraction', 
                     SCRNAFeatureExtractionStep(self.workdir, script_dir, config))
        self.add_step('tree_input', 
                     SCRNATreeInputStep(self.workdir, script_dir, config))
        self.add_step('tree_building', 
                     SCRNATreeBuildingStep(self.workdir, script_dir, config))
    
    def run(self, sample_id: str, mutation_list: Path, bam_file: Path,
            barcode_file: Path, metadata_file: Optional[Path] = None,
            read_len: int = 91, cellnum: int = 155, threads: int = 4,
            running_type: str = "benchmark", ase_filepath: Optional[str] = None,
            steps: List[str] = None, parallel: bool = False) -> Dict[str, Any]:
        """
        Run the complete scRNA pipeline
        
        Args:
            sample_id: Sample identifier
            mutation_list: Path to mutation list file
            bam_file: Path to BAM file
            barcode_file: Path to barcode file
            metadata_file: Path to metadata file (optional)
            read_len: Read length
            cellnum: Number of cells
            threads: Number of threads
            running_type: Running type (benchmark/production)
            ase_filepath: ASE file path (optional)
            steps: List of steps to run (None for all)
            parallel: Whether to run feature extraction and tree input in parallel
            
        Returns:
            Dictionary with all step results
        """
        self.logger.info(f"Starting scRNA pipeline for sample {sample_id}")
        
        # Prepare arguments for each step
        step_kwargs = {
            'feature_extraction': {
                'sample_id': sample_id,
                'mutation_list': mutation_list,
                'bam_file': bam_file,
                'barcode_file': barcode_file,
                'read_len': read_len,
                'running_type': running_type,
                'ase_filepath': ase_filepath,
                'threads': threads
            },
            'tree_input': {
                'sample_id': sample_id,
                'mutation_list': mutation_list,
                'bam_file': bam_file,
                'barcode_file': barcode_file,
                'cellnum': cellnum,
                'threads': min(threads, 2)  # Tree input may only need 2 threads
            },
            'tree_building': {
                'sample_id': sample_id,
                'cellnum': cellnum,
                'metadata_file': metadata_file,
                'barcode_file': barcode_file,
                # Features file from tree input step
                'features_file': self.workdir / '02_treeinput' / 'treeinput' / 'features_file.txt'
            }
        }
        
        # Determine steps to run
        steps_to_run = steps or ['feature_extraction', 'tree_input', 'tree_building']
        
        # Run first two steps in parallel if requested
        if parallel and 'feature_extraction' in steps_to_run and 'tree_input' in steps_to_run:
            self.logger.info("Running feature extraction and tree input in parallel")
            
            with ThreadPoolExecutor(max_workers=2) as executor:
                futures = {}
                
                if 'feature_extraction' in steps_to_run:
                    futures['feature_extraction'] = executor.submit(
                        self.run_step, 'feature_extraction', 
                        **step_kwargs['feature_extraction']
                    )
                    steps_to_run.remove('feature_extraction')
                
                if 'tree_input' in steps_to_run:
                    futures['tree_input'] = executor.submit(
                        self.run_step, 'tree_input',
                        **step_kwargs['tree_input']
                    )
                    steps_to_run.remove('tree_input')
                
                # Collect results from parallel steps
                for name, future in futures.items():
                    self.results[name] = future.result()
        
        # Run remaining steps sequentially
        for step_name in steps_to_run:
            if step_name in self.steps:
                self.run_step(step_name, **step_kwargs.get(step_name, {}))
        
        # Add summary information
        self.results['summary'] = {
            'sample_id': sample_id,
            'workdir': str(self.workdir),
            'steps_completed': list(self.results.keys()),
            'tree_file': self.get_tree_file()
        }
        
        return self.results
    
    def get_tree_file(self) -> Optional[Path]:
        """
        Get the final tree file path
        
        Returns:
            Path to tree file or None if not found
        """
        return self.get_step_output('tree_building', 'tree_file')