#!/usr/bin/env python3

def download_shell():
    import subprocess
    return_code = subprocess.call("wget -r -np -nd -R 'index.html*' -e robots=off https://www1.ncdc.noaa.gov/pub/data/igra/data/data-por", shell=True)  

def main():
    import os
    #
    # IGRA source directory
    #
    DATADIR=/tmp/data/igra
    
    # create directories
    if not os.path.isdir(DATADIR):
        os.makedirs(DATADIR)
    
    # Download everything
    status = download_shell()
    if status != 0:
        raise RuntimeError("Download interrupted")
    


if __name__ == "__main__":
    main()
