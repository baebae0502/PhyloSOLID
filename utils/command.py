"""
Command execution tools for calling external scripts
"""

import subprocess
import logging
import shlex
import os
from pathlib import Path
from typing import Union, List, Optional, Dict, Any

logger = logging.getLogger(__name__)

class CommandRunner:
    """Command executor for running external commands and scripts"""
    
    def __init__(self, workdir: Union[str, Path] = None, env: Dict[str, str] = None):
        """
        Initialize command runner
        
        Args:
            workdir: Working directory for commands
            env: Environment variables
        """
        self.workdir = Path(workdir) if workdir else Path.cwd()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.env = env or os.environ.copy()
        
    def run(self, cmd: Union[str, List[str]], check: bool = True, 
            capture_output: bool = False, **kwargs) -> subprocess.CompletedProcess:
        """
        Run a command
        
        Args:
            cmd: Command string or list of arguments
            check: Whether to check return code
            capture_output: Whether to capture stdout/stderr
            **kwargs: Additional arguments for subprocess.run
            
        Returns:
            subprocess.CompletedProcess object
        """
        if isinstance(cmd, list):
            cmd_str = ' '.join(shlex.quote(str(c)) for c in cmd)
        else:
            cmd_str = cmd
            
        logger.info(f"Running: {cmd_str[:200]}...")
        logger.debug(f"Working directory: {self.workdir}")
        
        try:
            result = subprocess.run(
                cmd_str,
                shell=True,
                cwd=self.workdir,
                env=self.env,
                capture_output=capture_output,
                text=True,
                check=check,
                **kwargs
            )
            
            if capture_output and result.stdout:
                logger.debug(f"STDOUT: {result.stdout[:500]}...")
            if capture_output and result.stderr:
                logger.debug(f"STDERR: {result.stderr[:500]}...")
                
            return result
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with exit code {e.returncode}")
            if e.stdout:
                logger.error(f"STDOUT: {e.stdout}")
            if e.stderr:
                logger.error(f"STDERR: {e.stderr}")
            raise

    def run_script(self, script_path: Union[str, Path], 
                   args: List[str] = None, **kwargs) -> subprocess.CompletedProcess:
        """
        Run a shell script
        
        Args:
            script_path: Path to script file
            args: List of arguments to pass to script
            **kwargs: Additional arguments for run()
            
        Returns:
            subprocess.CompletedProcess object
        """
        script_path = Path(script_path)
        if not script_path.exists():
            raise FileNotFoundError(f"Script not found: {script_path}")
            
        # Ensure script is executable
        script_path.chmod(0o755)
        
        cmd = [str(script_path)] + (args or [])
        return self.run(cmd, **kwargs)


class RScriptRunner(CommandRunner):
    """R script executor"""
    
    def __init__(self, rscript_path: str = "Rscript", **kwargs):
        """
        Initialize R script runner
        
        Args:
            rscript_path: Path to Rscript executable
            **kwargs: Arguments for CommandRunner
        """
        super().__init__(**kwargs)
        self.rscript = rscript_path
        
    def run(self, script_path: Union[str, Path], 
            args: List[str] = None, **kwargs) -> subprocess.CompletedProcess:
        """
        Run an R script
        
        Args:
            script_path: Path to R script
            args: Arguments to pass to script
            **kwargs: Additional arguments for run()
            
        Returns:
            subprocess.CompletedProcess object
        """
        script_path = Path(script_path)
        if not script_path.exists():
            raise FileNotFoundError(f"R script not found: {script_path}")
            
        cmd = [self.rscript, str(script_path)] + (args or [])
        return super().run(cmd, **kwargs)
