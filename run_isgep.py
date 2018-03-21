"""Place this script in any of the example folder to run ISGEP"""

import os

working_dir = os.getcwd()
configuration = os.path.join(os.getcwd(),"ui.conf")
isgep_path = os.path.join(os.path.dirname(os.path.dirname(working_dir)),"src/isgep_master.py")
os.system("python "+isgep_path+" "+configuration)
