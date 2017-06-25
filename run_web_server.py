# To launch the GUI:
# On ubuntu/linux, run:
# > cd /path/to/this/directory/
# > sudo screen python run_web_server.py -S "Pandexo Server"
# This will enable the web server to run as a background process.

import pandexo.engine.run_online as ro
ro.main()