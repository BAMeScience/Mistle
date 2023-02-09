#time ./mistle-build -i '/home/ynowatzk/data/benchmarking/9MM/msp/9MM_DB.msp /home/ynowatzk/data/benchmarking/CRAP/msp/crap.msp /home/ynowatzk/data/benchmarking/human/msp/human.msp' -o /home/ynowatzk/data/benchmarking/9MM/study/mistle/entrapment/ -n 64 -t 8
#time ./mistle-build -i '/home/ynowatzk/data/benchmarking/9MM/msp/9MM_DB_decoy.msp /home/ynowatzk/data/benchmarking/CRAP/msp/crap_decoy.msp /home/ynowatzk/data/benchmarking/human/msp/human_decoy.msp' -o /home/ynowatzk/data/benchmarking/9MM/study/mistle/entrapment_decoy/ -n 64 -t 8


time ./mistle-search -s ~/data/benchmarking/9MM/mgf/9MM_FASP.mgf -i ~/data/benchmarking/9MM/study/mistle/entrapment/ -o 9MM_FASP.pin -p 10.0 -b 0.2 -t 8
time ./mistle-search -s ~/data/benchmarking/9MM/mgf/9MM_PPID.mgf -i ~/data/benchmarking/9MM/study/mistle/entrapment/ -o 9MM_PPID.pin -p 10.0 -b 0.2 -t 8
time ./mistle-search -s ~/data/benchmarking/9MM/mgf/9MM_Run_1.mgf -i ~/data/benchmarking/9MM/study/mistle/entrapment/ -o 9MM_Run_1.pin -p 10.0 -b 0.2 -t 8
time ./mistle-search -s ~/data/benchmarking/9MM/mgf/9MM_Run_2.mgf -i ~/data/benchmarking/9MM/study/mistle/entrapment/ -o 9MM_Run_2.pin -p 10.0 -b 0.2 -t 8


time ./mistle-search -s ~/data/benchmarking/9MM/mgf/9MM_FASP.mgf -i ~/data/benchmarking/9MM/study/mistle/entrapment_decoy/ -o 9MM_FASP.decoy.pin -p 10.0 -b 0.2 -t 8
time ./mistle-search -s ~/data/benchmarking/9MM/mgf/9MM_PPID.mgf -i ~/data/benchmarking/9MM/study/mistle/entrapment_decoy/ -o 9MM_PPID.decoy.pin -p 10.0 -b 0.2 -t 8
time ./mistle-search -s ~/data/benchmarking/9MM/mgf/9MM_Run_1.mgf -i ~/data/benchmarking/9MM/study/mistle/entrapment_decoy/ -o 9MM_Run_1.decoy.pin -p 10.0 -b 0.2 -t 8
time ./mistle-search -s ~/data/benchmarking/9MM/mgf/9MM_Run_2.mgf -i ~/data/benchmarking/9MM/study/mistle/entrapment_decoy/ -o 9MM_Run_2.decoy.pin -p 10.0 -b 0.2 -t 8



python ../scripts/merge_pin_output.py -t ./9MM_FASP.pin -d ./9MM_FASP.decoy.pin -o ./9MM_FASP.td.pin --update_delta_scores
python ../scripts/merge_pin_output.py -t ./9MM_PPID.pin -d ./9MM_PPID.decoy.pin -o ./9MM_PPID.td.pin --update_delta_scores
python ../scripts/merge_pin_output.py -t ./9MM_Run_1.pin -d ./9MM_Run_1.decoy.pin -o ./9MM_Run_1.td.pin --update_delta_scores
python ../scripts/merge_pin_output.py -t ./9MM_Run_2.pin -d ./9MM_Run_2.decoy.pin -o ./9MM_Run_2.td.pin --update_delta_scores

