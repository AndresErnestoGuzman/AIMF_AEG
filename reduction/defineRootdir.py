#defineRootdir.py
## Define almaimf_rootdir and put it in path
if not 'almaimf_rootdir' in locals():
    if 'ALMAIMF_ROOTDIR' in os.environ:
        almaimf_rootdir = os.environ['ALMAIMF_ROOTDIR']
    elif not os.getenv('ALMAIMF_ROOTDIR') is None:
        almaimf_rootdir = os.environ['ALMAIMF_ROOTDIR'] = os.getenv('ALMAIMF_ROOTDIR')
    else:
        try:
            import metadata_tools
            almaimf_rootdir = os.environ['ALMAIMF_ROOTDIR'] = os.path.split(metadata_tools.__file__)[0]
        except ImportError:
            raise ValueError("metadata_tools not found on path; make sure to "
                             "specify ALMAIMF_ROOTDIR environment variable "
                             "or your PYTHONPATH variable to include the directory"
                             " containing the ALMAIMF code.")
else:
    os.environ['ALMAIMF_ROOTDIR'] = almaimf_rootdir
if not almaimf_rootdir in sys.path:
    sys.path.append(os.getenv('ALMAIMF_ROOTDIR'))
    