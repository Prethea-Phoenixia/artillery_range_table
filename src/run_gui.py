"""
This script is required to help PyInstaller to locate the exterior_ballistics package.
And pretty useful for debugging.
"""

from exterior_ballistics.gui import main
from multiprocessing import freeze_support


if __name__ == "__main__":
    freeze_support()  # required for Windows
    main()
