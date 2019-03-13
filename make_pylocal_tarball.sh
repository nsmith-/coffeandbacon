tar -czf pylocal.tgz -C ~/.local/lib/python3.6/ site-packages
xrdcp pylocal.tgz root://cmseos.fnal.gov//store/user/$USER/pylocal.tgz
