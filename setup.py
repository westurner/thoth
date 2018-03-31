import os
import sys

argc = len(sys.argv)

if (argc == 1):
   print('Unknown option: please specify if you want to install or uninstall the package')
   exit()
elif (argc == 2):
   option = str(sys.argv[1])
   option = option.rstrip()
   if (option == 'install'):
      print('WARNING!!! Install THOTH package into default python directories')
      command = 'python setup_cfg.py install'
   elif (option == 'uninstall'):
      try:
         with open('installation_path.save'): 
            pass
            f = open('installation_path.save','r')
            installation_path = f.readline()
            f.close()
            #print(installation_path)
            os.system('rm -rf ' + installation_path)
            os.system('rm -rf build/ installation_path.save')
            print('Uninstall complete.')
      except IOError:
         exception_error = -1;
         print("Cannot find 'installation_path.save'")
   else:
      print('ERROR!!! Unrecognized option')
      exit()
elif (argc == 3):
   option = str(sys.argv[1])
   option = option.rstrip()
   if (option == 'install'):
      installation_path = str(sys.argv[2])
      installation_path = installation_path.rstrip()
      print(installation_path)
      print('Installation path: ' + installation_path)
      command = 'python setup_cfg.py install --prefix=' + installation_path + ' --install-lib=' + installation_path + ' --install-platlib=' + installation_path + ' --install-scripts=' + installation_path + ' --install-data=' + installation_path
      f = open('installation_path.save','w')
      f.write(installation_path)
      f.flush()
      f.close()
      #print(command)
      os.system(command)
   elif (option == 'uninstall'):
      print('WARNING!!! Argument not used.')
      f = open('installation_path.save','r')
      installation_path = f.readline()
      f.close()
      #print(installation_path)
      os.system('rm -rf ' + installation_path)
      os.system('rm -rf build/ installation_path.save')
      print('Uninstall complete.')
   else:
      print('ERROR!!! Unrecognized option')
      exit()
