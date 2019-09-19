from parsl.providers import SlurmProvider
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.addresses import address_by_hostname

x509_proxy = 'x509up_u%s'%(os.getuid())

wrk_init = '''
export XRD_RUNFORKHANDLER=1
export X509_USER_PROXY=$WORK/%s
'''%(x509_proxy)

twoGB = 2048
nproc = 16

sched_opts = '''
#SBATCH --cpus-per-task=%d
#SBATCH --mem-per-cpu=%d
''' % (nproc, twoGB, )


slurm_htex = Config(
    executors=[
        HighThroughputExecutor(
            label="coffea_parsl_slurm",
            address=address_by_hostname(),
            prefetch_capacity=0,
            max_workers=nproc,
            provider=SlurmProvider(
                launcher=SrunLauncher(),
                init_blocks=4,
                max_blocks=4,
                nodes_per_block=1,
                partition='batch,guest,gpu',
                scheduler_options=sched_opts,   # Enter scheduler_options if needed
                worker_init=wrk_init,         # Enter worker_init if needed
                walltime='02:00:00'
            ),
        )
    ],
    strategy=None,
)
