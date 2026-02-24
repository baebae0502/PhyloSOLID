"""
File and directory utility functions
"""

import os
import shutil
from pathlib import Path
from typing import List, Union, Optional

def ensure_dir(path: Union[str, Path]) -> Path:
    """
    Ensure directory exists, create if necessary
    
    Args:
        path: Directory path
        
    Returns:
        Path object
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path

def find_files(directory: Union[str, Path], pattern: str) -> List[Path]:
    """
    Recursively find files matching pattern
    
    Args:
        directory: Directory to search
        pattern: File pattern (e.g., "*.txt")
        
    Returns:
        List of matching file paths
    """
    directory = Path(directory)
    return list(directory.rglob(pattern))

def copy_file(src: Union[str, Path], dst: Union[str, Path]) -> Path:
    """
    Copy file from source to destination
    
    Args:
        src: Source file path
        dst: Destination file path
        
    Returns:
        Destination path
    """
    src = Path(src)
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst
