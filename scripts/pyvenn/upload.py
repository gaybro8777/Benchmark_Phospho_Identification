from imp import load_source
from subprocess import call

meta = load_source("_", "./setup.py").meta
tarball = "dist/{}-{}.tar.gz".format(meta.__name__, meta.__version__)

call(["python", "setup.py", "sdist"])
call(["twine", "upload", "-r", "pypi", tarball])
