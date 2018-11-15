# Hall A CompMon

## Compiling
In the directory src/config.d createa  file like so:
```bash
BINDIR=<path_to_output_executable>
EVIO_LIBPATH=<basedir_of_local_evio_env>                                    
EVIO_INCLUDEPATH=\${EVIO}/include                                               
EVIO_LIB=lib/libevio.so
```
Then in in the src/ directory:
```bash
./config config.d/<your_config_file>
make
```

## Replaying raw data files
By default compmon looks for files in /data/cmu, where the raw data files are prefixed with FadcCalo2016. You can change this by creating a copy of the 'compmon.config' file in the base directory of this repository. You can change the following variables:
```
  dataPath 		# Path to raw data files
  dataPrefix	# Prefix of raw data file, in case they don't start with FadcCalo2016_
  outPath 		# Path to where ROOT files will be stored
  outPrefix		# Rootfiles will be prefixed with this string
```

To replay a given run (i.e. 2899), (assuming the compmon binary is in your path):
```
   compmon -r 2899 -c <your config file>
```
Optionally, to limit the number of entries to replay you can use ``-n <num_entries>``. 
