system("rm -rf /mnt/scratch_2/output/230124_VL00165_23_AACHJT5M5_mapping/out/consensus_fa1/");
system("mkdir /mnt/scratch_2/output/230124_VL00165_23_AACHJT5M5_mapping/out/consensus_fa1/");
system("ls /mnt/scratch_2/output/230124_VL00165_23_AACHJT5M5_mapping/out/consensus_fa/*.gz > consensus_fasta_list.txt");


$file = "consensus_fasta_list.txt";

open(in, $file);

while (<in>){
    chomp;
    $s = $_;
    @simple = split/consensus\_fa\//,$s;
    for $m (0..$#simple){}
    #print "$simple[0]\n";
    if ($simple[0] !~ "neg") {
        #print "$simple[0]\t$simple[1]\n";
        @label = split/\_/,$simple[1];
        for $l (0..$#label){}
        #print "${s}++\t$label[0]\n";
        #print "$label[0]\n";
        system("cp ${s} /mnt/scratch_2/output/230124_VL00165_23_AACHJT5M5_mapping/out/consensus_fa1/$label[0].fa.gz");
        if($label[0] !~ "neg"){
        print "$label[0]\t/mnt/scratch_2/output/230124_VL00165_23_AACHJT5M5_mapping/out/consensus_fa1/$label[0].fa.gz\n";
        }      
    }   
}               

close(in);

exit;