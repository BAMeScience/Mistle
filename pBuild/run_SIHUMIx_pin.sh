#This one##time -p ./mistle-build -i '/home/ynowatzk/data/benchmarking/SIHUMIx/msp/SIHUMIx.msp /home/ynowatzk/data/benchmarking/human/msp/human.msp' -o /home/ynowatzk/data/benchmarking/SIHUMIx/study/mistle/entrapment/ -n 64 -t 8
#time ./mistle-build -i '/home/ynowatzk/data/benchmarking/SIHUMIx/msp/SIHUMIx.msp' -o /home/ynowatzk/data/benchmarking/SIHUMIx/study/mistle/target/ -n 64 -t 8
#THIS ONE##time -p ./mistle-build -i '/home/ynowatzk/data/benchmarking/SIHUMIx/msp/SIHUMIx_decoy.msp /home/ynowatzk/data/benchmarking/human/msp/human_decoy.msp' -o /home/ynowatzk/data/benchmarking/SIHUMIx/study/mistle/entrapment_decoy/ -n 64 -t 8
#time ./mistle-build -i '/home/ynowatzk/data/benchmarking/SIHUMIx/msp/SIHUMIx_decoy.msp' -o /home/ynowatzk/data/benchmarking/SIHUMIx/study/mistle/decoy/ -n 64 -t 8

#time ./mistle-search -s ~/data/benchmarking/SIHUMIx/mgf/S05.mgf -i ~/data/benchmarking/SIHUMIx/study/mistle/target/ -o S05.csv -p 10.0 -b 0.02 -t 8
#time ./mistle-search -s ~/data/benchmarking/SIHUMIx/mgf/S05.mgf -i ~/data/benchmarking/SIHUMIx/study/mistle/decoy/ -o S05.csv -p 10.0 -b 0.02 -t 8
#time ./mistle-search -s ~/data/benchmarking/SIHUMIx/mgf/S06.mgf -i ~/data/benchmarking/SIHUMIx/study/mistle/decoy/ -o S06.csv -p 10.0 -b 0.02 -t 8
#time ./mistle-search -s ~/data/benchmarking/SIHUMIx/mgf/S06.mgf -i ~/data/benchmarking/SIHUMIx/study/mistle/target/ -o S06.csv -p 10.0 -b 0.02 -t 8


time ./mistle-search -s ~/data/benchmarking/SIHUMIx/mgf/S05.mgf -i ~/data/benchmarking/SIHUMIx/study/mistle/entrapment_decoy/ -o S05.csv -p 10.0 -b 0.02 -t 8
time ./mistle-search -s ~/data/benchmarking/SIHUMIx/mgf/S06.mgf -i ~/data/benchmarking/SIHUMIx/study/mistle/entrapment_decoy/ -o S06.csv -p 10.0 -b 0.02 -t 8
time ./mistle-search -s ~/data/benchmarking/SIHUMIx/mgf/S05.mgf -i ~/data/benchmarking/SIHUMIx/study/mistle/entrapment/ -o S05.csv -p 10.0 -b 0.02 -t 8
time ./mistle-search -s ~/data/benchmarking/SIHUMIx/mgf/S06.mgf -i ~/data/benchmarking/SIHUMIx/study/mistle/entrapment/ -o S06.csv -p 10.0 -b 0.02 -t 8


