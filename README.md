# Using UprootJob in LCG94python3
On any machine that has access to `/cvmfs`, the following setup script will install the relevant packages as user.
```bash
source setup_lcg.sh
```
On future use, run `source env_lcg.sh`, and start a notebook server with, e.g. `jupyter notebook --no-browser --port $PORT`, substituing your favorite port: pick a random integer in (8000,65535).
Often, you'll need an ssh tunnel, which can be accomplished via, e.g. `ssh -L $PORT:localhost:$PORT server.address`

# Using a docker container
The docker container can run both uproot and striped jobs (although xrootd has not been tested yet.)
To run, execute
```bash
./run_container.sh
```
then go to http://localhost:8888/ and insert the token
