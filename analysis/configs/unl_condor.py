from parsl.providers import CondorProvider
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_hostname

x509_proxy = 'x509up_u%s'%(os.getuid())

wrk_init = '''
source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-slc6-gcc62-opt/setup.sh
export PATH=`pwd`/.local/bin:$PATH
export PYTHONPATH=`pwd`/.local/lib/python3.6/site-packages:$PYTHONPATH

export X509_USER_PROXY=`pwd`/%s
mkdir -p ./coffea_parsl_condor
'''%(x509_proxy)

twoGB = 1024
nproc = 16

condor_cfg = '''
universe = Docker

+WantWholeNode = True
use_x509userproxy = true
+WantDocker = True
docker_image = "opensciencegrid/osgvo-el6"

transfer_output_files = coffea_parsl_condor
''' #% (nproc, ) # twoGB*nproc, 

#RequestMemory = %d
#RequestCpus = %d
#RequestDisk = 1048576

xfer_files = ['%s/.local' % (os.environ['HOME'], ), '/tmp/%s' % (x509_proxy, )]

config = Config(
    executors=[
        HighThroughputExecutor(
            label="coffea_parsl_condor",
            address=address_by_hostname(),
            prefetch_capacity=0,
            cores_per_worker=1,
            max_workers=nproc,
            worker_logdir_root='./',
            provider=CondorProvider(
                init_blocks=8,
                max_blocks=200,
                nodes_per_block=1,
                worker_init = wrk_init,
                transfer_input_files=xfer_files,
                scheduler_options=condor_cfg
            ),
        )
    ],
    strategy=None,
)

