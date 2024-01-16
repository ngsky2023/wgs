#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#---------------------------------------------------------
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
my $ver="2.0";
my (%tr_hgvsc,%tr_hgvsp,%tr_hgvsc_noRS,%tr_hgvsp_noRS,%gene_chain,%hs_nm,%left_align,%hs_idinfo,%hs_idinfo_noRS,%hs_NM,%gene_mainNM,%hs_idinfo_hgvsc,%hs_idinfo_noRS_hgvsc,);

my $main_NM = "$Bin/gene_NM.list";
open NM,"<$main_NM" or die $!;
while (<NM>){
	chomp;
	next unless /\t/;
	my @tmp = split /\t/,$_;
	$gene_mainNM{$tmp[0]} = $tmp[1];
}
close NM;
open TR,"<$ARGV[1]" or die $!;
my $n="0";
while (<TR>){
    chomp;
    next if $_ =~ /^#/;
    my @tr = split /\t/, $_;
    my ($chr,$pos,$rs,$ref,$alt,) = ($tr[0],$tr[1],$tr[2],$tr[3],$tr[4]);
  
    my ($NM,$gene,$chain,$variant_info,$left_align_info,$exinfo,)=($tr[11],$tr[12],$tr[13],$tr[14],$tr[16],$tr[15],);
    if ($NM =~ /(.M.*?)\s+\(.*?\)/){
        $NM = $1;
    }
    my $hgvsc = "";
    my $hgvsp = "";
	my $exid = "";
    my $left_align_hgvsc;
    if ($variant_info =~ /\/(c\..*?)\/(.*)/){   #chr10:g.43610119G>A/c.2071G>A/p.G691S
        $hgvsc = $1;
		 $hgvsp = $2;
#		if ($hgvsc =~ /^c.(\d+)([ACGT])>([ACGTYRN])$/){
#			$hgvsc = "c\.$2$1$3";
#		}
        if($hgvsp =~ /(.*?)\*fs\*\d+$/){
            $hgvsp = $1. "Xfs";
		}elsif ($hgvsp =~ /(.*)\*\d+$/){
            $hgvsp = $1;
        }elsif($hgvsp =~ /(.*)\*$/){
            $hgvsp = $1. "X";
        }
    }
    if ($left_align_info =~ /\;left_align_cDNA=(c\..*?);/){
        $left_align_hgvsc = $1;    
    }else{$left_align_hgvsc = $hgvsc}
    if ($exinfo =~ /inside\_\[(.*)\]/){
		$exid = $1;
	}
    my @alts = split /,/, $alt;
    for my $key (@alts){
        my $delsize = length($ref) - length($key);
		if (length($ref) > length($key) && (substr($ref, 0, length($ref)-$delsize) eq $key) ) {
            $pos = $pos + (length($ref)-$delsize);
            $ref = substr($ref, length($ref)-$delsize, $delsize);
            $alt = "-" ;
        }elsif (length($ref) > length($key) && (substr($ref, 0, 1 ) eq substr($key, 0, 1)) ){
            $pos = $pos + 1;
            $ref = substr($ref, 1, length($ref)-1);
            $alt = substr($key, 1, length($key)-1); 
        }elsif ((length($ref) >1 && length($key)>1) && length($ref) == length($key) && (substr($ref, -1, length($key)-1)  eq substr($key,-1,length($ref) -1)) ) {
            $pos = $pos;
            $ref = substr($ref, 0, length($ref)-1);
            $alt = substr($key, 0 ,length($key)-1) ;
        }elsif ((length($ref) >1 && length($key)>1) && length($ref) == length($key) && (substr($ref, 0 , 1)  eq  substr($key, 0, 1)) && (substr($ref, -1, 1)  ne  substr($key, -1, 1)) ) {
            $pos = $pos + 1;
            $ref = substr($ref, 1 ,length($ref)-1);
            $alt = substr($key, 1 ,length($key)-1) ;
        }elsif (length($ref) < length($key) && (substr($key, 0, length($key)+$delsize) eq $ref) ) {
		#$pos = $pos + (length($alt)+$delsize);
		$alt = substr($key, length($key)+$delsize, (-$delsize));
		$ref = "-" ;
        }elsif (length($ref) < length($key) && ( substr($ref, 0, length($key)-length($ref)) eq substr($key, 0, length($key)-length($ref)) ) ) {
		$pos = $pos - $delsize;
		$alt = substr($key, (-$delsize), length($key)-$delsize);
		$ref = substr($ref, (-$delsize), length($ref)-$delsize);
        }elsif (length($ref) < length($key) && ( substr($ref, 0, 1) eq substr($key, 0 , 1) ) ) {   #yangrutao 20180206
                $pos = $pos + 1;
                $alt = substr($key, 1, length($key)-1);
                $ref = substr($ref, 1, length($ref)-1);
        }elsif (length($ref) == 1 && length($key) == 1 ) {
            $pos = $pos;
            $alt = $alt;
        }
		
		if ($rs ne ".") {
			$tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$NM}{$rs} = $hgvsc;
			$tr_hgvsp{$chr}{$pos}{$ref}{$alt}{$NM}{$rs} = $hgvsp;
		}
		
		$tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$NM} = $hgvsc;
		$tr_hgvsp_noRS{$chr}{$pos}{$ref}{$alt}{$NM} = $hgvsp;
		$left_align{$chr}{$pos}{$ref}{$alt}{$NM}{$rs} = $left_align_hgvsc;
		$gene_chain{$gene} = $chain;
		$hs_nm{$gene}{$NM} = $NM;
		
		$hs_idinfo{$chr}{$pos}{$ref}{$alt}{$NM}{$rs}="$gene:$NM:$exid:$hgvsc:$hgvsp";
		$hs_idinfo_noRS{$chr}{$pos}{$ref}{$alt}{$NM}="$gene:$NM:$exid:$hgvsc:$hgvsp";
		$hs_idinfo_hgvsc{$chr}{$pos}{$ref}{$alt}{$NM}{$rs}="$gene:$NM:$exid:$hgvsc";
		if ($exid =~ /intron_between_exon_(\d+)_and_(\d+)/){
			my ($exid1 , $exid2) = ($1 , $2);
			if ($hgvsc =~ /\+/){
				$exid = "exon$exid1";
			}else{
				$exid = "exon$exid2";
			}
		}
		$hs_idinfo_noRS_hgvsc{$chr}{$pos}{$ref}{$alt}{$NM}="$gene:$NM:$exid:$hgvsc";
		
		$n++;
		$hs_NM{$NM}=$gene;
	}
	#if ($pos eq "94212931"){print "$chr\t$pos\t$ref\t$alt\t$NM\t$rs\n";}
}
close TR;

open OUT,">$ARGV[2]" or die $!;
open LOG,">log" or die $!;
open IN,"<$ARGV[0]" or die $!;
while (<IN>){
	chomp;
    if ($_ =~ /^SampleName/){
        print OUT "$_\n";
        next;
    };
	my @dat=split /\t/,$_;
    my ($chr,$pos,$ref,$alt,) = ($dat[1],$dat[2],$dat[4],$dat[5]);
	my $rs = $dat[-1];
	#if ($dat[23] =~ /^rs/){
	#	$rs = $dat[23]; 
	#}else{
	#	print LOG "$chr\t$pos\t$ref\t$alt\t--\tNo RS\*!!\n";
	#}
    my $gene = $dat[12];
	my @gcc = split /[,;]/ , $gene;
	if ($#gcc > 0){
		$gene = (sort {length($a)<=>length($b)} @gcc)[0];
	}
	my $id_info;
    if ($gene_chain{$gene} eq "+"){
        my @info;	
        if ($dat[11] eq "exonic" || $dat[11] eq "exonic;exonic" || $dat[11] eq "exonic;splicing" ){
            my @var1 = split/,/, $dat[15];
            for (my $i=0;$i<15;$i++){
                print OUT "$dat[$i]\t";
            }
			my $info_NM;
            foreach my $key_var (@var1){
                @info = split /\:/,$key_var;
				if($info[1] =~ /(NM_\d+)\.\d/){
					$info_NM = $1;
				}else{
					$info_NM = $info[1];
				}
                if (exists $hs_nm{$gene}{$info_NM} and $info_NM eq $hs_nm{$gene}{$info_NM} && defined $rs  && defined  $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs} && defined $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs}){
                    print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t$rs\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
                    $info[3] = $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs};
                    $info[4] = $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs};
                    print OUT join ("\:",@info) ;
					print OUT ",";
                    print LOG "$info[3]\t$info[4]\n";
                }elsif (exists $hs_nm{$gene}{$info_NM} and $info_NM eq $hs_nm{$gene}{$info_NM} ){  #&& !defined $rs
                    print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t.\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
                    $info[3] = $tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$info_NM};
                    $info[4] = $tr_hgvsp_noRS{$chr}{$pos}{$ref}{$alt}{$info_NM};
                    print OUT join ("\:",@info) ;
					print OUT ",";
                    print LOG "$info[3]\t$info[4]\n";
                }else{
                    print OUT join ("\:",@info);
					print OUT ",";
                }
                
            }
			print OUT "\t";
            for (my $i=16;$i<$#dat;$i++){
                print OUT "$dat[$i]\t";
            }
            print OUT "$dat[-1]\n";
            
        }elsif($dat[11] eq "splicing" ){    #&& $hs_nm{$gene} =~ /$dat[13]/   -->test
			if ($dat[13] ne "."){
            #NM_000027:exon1:c.127+1->ATGCGG,NM_001171988:exon1:c.127+1->ATGCGG
            #NM_001195302:exon4:c.450-2A>C,NM_024598:exon5:c.504-2A>C
			my @info;
            for (my $i=0;$i<13;$i++){
                    print OUT "$dat[$i]\t";
            }	
				my $info_NM;
				my @var2 = split/\,/,$dat[13];
				foreach my $key_var (@var2){
				    @info = split /\:/,$key_var;
					if($info[0] =~ /(NM_\d+)\.\d/){
						$info_NM = $1;
					}else{
						$info_NM = $info[0];
					}
				    if ($info_NM eq $hs_nm{$gene}{$info_NM} && defined $rs && defined  $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs}){
				        print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t$rs\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
				        $info[2] = $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs};
				        #$info[4] = $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{info_NM};
				        print OUT join ("\:",@info);
						print OUT ",";
				        print LOG "$info[2]\t\-\t\n";
				    }elsif ($info_NM eq $hs_nm{$gene}{$info_NM} ){ #&& !defined $rs
				        print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t.\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
				        $info[2] = $tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$info_NM};
				        #$info[4] = $tr_hgvsp_noRS{$chr}{$pos}{$ref}{$alt}{info_NM};
				        print OUT join ("\:",@info);
						print OUT ",";
				        print LOG "$info[2]\t\-\t\n";
				    }else{
				        print OUT join ("\:",@info);
						print OUT ",";
				    }
				}
				print OUT "\t$dat[14]\t";
				foreach my $key_var (@var2){
				    @info = split /\:/,$key_var;
					if($info[0] =~ /(NM_\d+)\.\d/){
						$info_NM = $1;
					}else{
						$info_NM = $info[0];
					}
				    if ($info_NM eq $hs_nm{$gene}{$info_NM} && defined $rs && defined  $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs}){
				        print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t$rs\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
				        $info[2] = $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs};
				        #$info[4] = $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{info_NM};
				        print OUT join ("\:",@info);
						print OUT ",";
				        print LOG "$info[2]\t\-\t\n";
				    }elsif ($info_NM eq $hs_nm{$gene}{$info_NM} ){ #&& !defined $rs
				        print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t.\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
				        $info[2] = $tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$info_NM};
				        #$info[4] = $tr_hgvsp_noRS{$chr}{$pos}{$ref}{$alt}{info_NM};
				        print OUT join ("\:",@info);
						print OUT ",";
				        print LOG "$info[2]\t\-\t\n";
				    }else{
				        print OUT join ("\:",@info);
						print OUT ",";
				    }
				}
			}else{
				my @info;
				for (my $i=0;$i<13;$i++){
                    print OUT "$dat[$i]\t";
				}
				print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t$rs\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
				my $info_tmp = $hs_idinfo_noRS_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
				#$info[1] = $gene_mainNM{$geng};
				my $info_tr = $tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
				#$info[4] = $tr_hgvsp_noRS{$chr}{$pos}{$ref}{$alt}{info_NM};
				#print OUT join ("\:",@info);
				print OUT "$info_tmp\t\.\t$info_tmp,";
				#print OUT ",";
				#print LOG "$info[2]\t\-\t\n";
				print LOG "$info_tr\t-\t\n";
			}
			print OUT "\t";
            for (my $i=16;$i<$#dat;$i++){
                    print OUT "$dat[$i]\t";
            }
             print OUT "$dat[-1]\n";
            
        }elsif ($dat[11] eq "intronic"){
			my $id_info="";
            for (my $i=0;$i<15;$i++){
                print OUT "$dat[$i]\t";
            }
			
            #if ($dat[15] eq "."){
					if ( defined $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs}){
						if ($tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs} =~ /c\.\d+(\+|\-)(\d+)[A-Z]\>[A-Z]/){
							if ($2 <=10){
								$id_info = $hs_idinfo{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs};
							}
						}elsif($tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs} =~ /c\.\d+(\+|\-)(\d+)_\d+(\+|\-)(\d+)(del|ins|dup)/){ #c.4001+2_4001+5delTAAC
							if ($2 <=10){
								$id_info = $hs_idinfo{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs};
							}
						}elsif($tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs} =~ /c\.\d+(\+|\-)(\d+)(del|ins|dup)/){
							if ($2 <=10){
								$id_info = $hs_idinfo{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs};
							}
						}
						print OUT $id_info;
						print OUT ",";
					}elsif (defined $tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}} ){ 
						if ($tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}} =~ /c\.\d+(\+|\-)(\d+)[A-Z]\>[A-Z]/){
							if ($2 <=10){
								$id_info = $hs_idinfo_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
							}
						}elsif($tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}} =~ /c\.\d+(\+|\-)(\d+)_\d+(\+|\-)(\d+)(del|ins|dup)/){  
							if ($2 <=10){
								$id_info = $hs_idinfo_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
							}
						}elsif($tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}} =~ /c\.\d+(\+|\-)(\d+)(del|ins|dup)/){  #c.2865-4delT
							if ($2 <=10){
								$id_info = $hs_idinfo_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
							}
						}
						print OUT $id_info;
						print OUT ",";
					}else{
						print OUT $dat[15];
					}
			print OUT "\t";
            for (my $i=16;$i<$#dat;$i++){
                print OUT "$dat[$i]\t";
            }
            print OUT "$dat[-1]\n";
            
        }else{
            print  OUT "$_\n";
        }
    }elsif ($gene_chain{$gene} eq "-"){
        my @info;
        if ($dat[11] eq "exonic" || $dat[11] eq "exonic;exonic" || $dat[11] eq "exonic;splicing" ) {  #&& $dat[15] =~ /$hs_nm{$gene}{$i}/
			my @all_nm;
			my $info_NM;
            my @var1 =split/\,/,$dat[15];
            for (my $i=0;$i<15;$i++){
                print OUT "$dat[$i]\t";
            }
            foreach my $key_var (@var1){
                @info = split /\:/,$key_var;
				if($info[1] =~ /(NM_\d+)\.\d/){
					$info_NM = $1;
				}else{
					$info_NM = $info[1];
				}
                if (exists  $hs_nm{$gene}{$info_NM} and $info_NM eq  $hs_nm{$gene}{$info_NM} && defined $rs && defined $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs} && defined $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs}){
                    print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t$rs\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
                    $info[3] = $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs};
                    $info[4] = $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{$info_NM}{$rs};
                    print OUT join ("\:",@info) ;
					print OUT ",";
                    print LOG "$info[3]\t$info[4]\n";
                }elsif (exists $hs_nm{$gene}{$info_NM} and $info_NM eq  $hs_nm{$gene}{$info_NM}){ # && !defined $rs
                    print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t.\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment:\t";
                    $info[3] = $tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$info_NM};
                    $info[4] = $tr_hgvsp_noRS{$chr}{$pos}{$ref}{$alt}{$info_NM};
                    print OUT join ("\:",@info) ;
					print OUT ",";
                    print LOG "$info[3]\t$info[4]\n";
                }else{
                    print OUT join ("\:",@info);
					print OUT ",";
                }
			}
			print OUT "\t";
            for (my $i=16;$i<$#dat;$i++){
                print OUT "$dat[$i]\t";
            }
            print OUT "$dat[-1]\n";
            
        }elsif($dat[11] eq "splicing" ){ #&& $hs_nm{$gene} =~ /$dat[13]/
            if ($dat[13] ne "."){
			my @info;
			for (my $i=0;$i<13;$i++){
                    print OUT "$dat[$i]\t";
            }
				my $info_NM;
				my @var2 = split/\,/,$dat[13];
				foreach my $key_var (@var2){
				    @info = split /\:/,$key_var;
					if($info[0] =~ /(NM_\d+)\.\d/){
						$info_NM = $1;
					}else{
						$info_NM = $info[0]
					}
						
				    if ($info_NM eq $hs_nm{$gene}{$info_NM}){
				        print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t$rs\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment\t";
				        #$info[2] = $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM};
				        #$info[4] = $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{$info_NM};
				        print OUT join ("\:",@info);
						print OUT ",";
				        print LOG "$info[2]\t\-\t\n";
				    }else{
				        print OUT join ("\:",@info);
						print OUT ",";
				    }
				  
				}
				print OUT "\t$dat[14]\t";
				foreach my $key_var (@var2){
				    @info = split /\:/,$key_var;
					if($info[0] =~ /(NM_\d+)\.\d/){
						$info_NM = $1;
					}else{
						$info_NM = $info[0]
					}
						
				    if ($info_NM eq $hs_nm{$gene}{$info_NM}){
				        print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t$rs\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment\t";
				        #$info[2] = $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM};
				        #$info[4] = $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{$info_NM};
				        print OUT join ("\:",@info);
						print OUT ",";
				        print LOG "$info[2]\t\-\t\n";
				    }else{
				        print OUT join ("\:",@info);
						print OUT ",";
				    }
				  
				}
			}else{
				my @info;
				for (my $i=0;$i<13;$i++){
                    print OUT "$dat[$i]\t";
				}
				print LOG "$dat[0]\t$gene\t$gene_chain{$gene}\t$dat[11]\t$info[1]\t$rs\tleft_alignment:\t$info[3]\t$info[4]\t--\tright_alignment\t";
				my $info_tr = $tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
				my $info_tmp = $hs_idinfo_noRS_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
				#$info[2] = $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$info_NM};
				#$info[4] = $tr_hgvsp{$chr}{$pos}{$ref}{$alt}{$info_NM};
				#print OUT join ("\:",@info);
				print OUT "$info_tmp\t\.\t$info_tmp,";
				#print OUT ",";
				#print LOG "$info[2]\t\-\t\n";
				print LOG "$info_tr\t\-\t\n";
			}
			print OUT "\t";
            for (my $i=16;$i<$#dat;$i++){
                    print OUT "$dat[$i]\t";
            }
             print OUT "$dat[-1]\n";
            
        }elsif ($dat[11] eq "intronic"){
			my $id_info="";
            for (my $i=0;$i<15;$i++){
                print OUT "$dat[$i]\t";
            }
			
            #if ($dat[15] eq "."){
					if (defined $tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs} ){
						if ($tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs} =~ /c\.\d+(\+|\-)(\d+)[A-Z]\>[A-Z]/){
							if ($2 <=10){
								$id_info = $hs_idinfo{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs};
							}
						}elsif($tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs} =~ /c\.\d+(\+|\-)(\d+)_\d+(\+|\-)(\d+)(del|ins|dup)/){ #c.4001+2_4001+5delTAAC
							if ($2 <=10){
								$id_info = $hs_idinfo{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs};
							}
						}elsif($tr_hgvsc{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs} =~ /c\.\d+(\+|\-)(\d+)(del|ins|dup)/){
							if ($2 <=10){
								$id_info = $hs_idinfo{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}}{$rs};
							}
						}
						print OUT $id_info;
						print OUT ",";
					}elsif (defined $tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}} ){ #&& !defined $rs
						if ($tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}} =~ /c\.\d+(\+|\-)(\d+)[A-Z]\>[A-Z]/){
							if ($2 <=10){
								$id_info = $hs_idinfo_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
							}
						}elsif($tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}} =~ /c\.\d+(\+|\-)(\d+)_\d+(\+|\-)(\d+)(del|ins|dup)/){  
							if ($2 <=10){
								$id_info = $hs_idinfo_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
							}
						}elsif($tr_hgvsc_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}} =~ /c\.\d+(\+|\-)(\d+)(del|ins|dup)/){  #c.2865-4delT
							if ($2 <=10){
								$id_info = $hs_idinfo_noRS{$chr}{$pos}{$ref}{$alt}{$gene_mainNM{$gene}};
							}
						}
						print OUT $id_info;
						print OUT ",";
					}else{
						print OUT $dat[15];
					}
			print OUT "\t";
            for (my $i=16;$i<$#dat;$i++){
                print OUT "$dat[$i]\t";
            }
            print OUT "$dat[-1]\n";
            
        }else{
            print  OUT "$_\n";
        }
    }else{
        print  OUT "$_\n";
    }
}

close IN;

close OUT;
close LOG;

my $DONE_TIME=time();
my $WORKDONE_TIME = &sub_format_datetime(localtime($DONE_TIME));
#print "Program Done Time:$WORKDONE_TIME\n";

#-------------------------------
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}



sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


