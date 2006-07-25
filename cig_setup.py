
import os, os.path, site, sys
from distutils.sysconfig import get_python_lib


cigprefix = os.path.expanduser(os.path.join("~", ".cig"))
cigpythondir = get_python_lib(plat_specific=False, prefix=cigprefix)
cigpyexecdir = get_python_lib(plat_specific=True, prefix=cigprefix)

check_pth_processing = None
write_script = None

sys_executable = os.path.normpath(sys.executable)
interpreter = sys_executable


def cig_check_pth_processing(self):
    if self.install_dir.startswith(cigprefix):
        return True
    return check_pth_processing(self)


def cig_write_script(self, script_name, contents, mode="t", blockers=()):
    if mode == "t":
        insertSitePatch = True
        newContents = []
        for line in contents.splitlines():
            if interpreter and line.startswith("#!"):
                line = "#!" + interpreter
            elif insertSitePatch and line.find('import ') != -1:
                newContents.append("import site")
                newContents.append("site.addsitedir('%s')" % cigpythondir)
                if cigpyexecdir != cigpythondir:
                    newContents.append("site.addsitedir('%s')" % cigpyexecdir)
                insertSitePatch = False
            newContents.append(line)
        contents = '\n'.join(newContents)
    return write_script(self, script_name, contents, mode, blockers)


def setup(**kwds):

    global interpreter
    interpreter = kwds.get('interpreter')
    if interpreter:
        del kwds['interpreter']

    # Create the CIG 'site-packages' directories, as needed.
    if not os.path.isdir(cigpythondir):
        os.makedirs(cigpythondir)
    if not os.path.isdir(cigpyexecdir):
        os.makedirs(cigpyexecdir)
    
    # Make '.pth' files work in the CIG 'site-packages' directories.
    site.addsitedir(cigpythondir)
    if cigpyexecdir != cigpythondir:
        site.addsitedir(cigpyexecdir)

    from ez_setup import use_setuptools
    use_setuptools()

    import setuptools

    # Patch setuptools (ugh).
    from setuptools.command.easy_install import easy_install
    global check_pth_processing, write_script
    check_pth_processing = easy_install.check_pth_processing
    write_script = easy_install.write_script
    easy_install.check_pth_processing = cig_check_pth_processing
    easy_install.write_script = cig_write_script

    # The "develop" command -- unlike "install" -- does not install
    # setuptools itself.  Work-around this problem.
    if setuptools.bootstrap_install_from:
        egg = setuptools.bootstrap_install_from
        setuptools.bootstrap_install_from = None
        setuptools.setup(script_args=['easy_install', egg])

    setuptools.setup(**kwds)
