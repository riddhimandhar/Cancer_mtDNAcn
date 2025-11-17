#!/usr/bin/perl -w 
#
## Process PCAWG analysis data and mtCopy number analysis data 
#

print "Reading Mito data ...\n"; 

open(FP,"../datasets/Cancer_mtGenome_analysis_NatGenet2020/TCMA-CopyNumber.tsv") or die; 
while($fp=<FP>)
{
    chomp($fp); 
    if($fp=~/^sample_id/) { next; }
    @arr=split(/\s+/,$fp);
    $listmtcn{$arr[0]}=$arr[2]; 
    undef(@arr);
}
close(FP); 

open(FP,"../datasets/Cancer_mtGenome_analysis_NatGenet2020/TCMA-MutationSNV.tsv") or die; 
while($fp=<FP>)
{
    chomp($fp); 
    if($fp=~/^sample_id/) { next; }
    @arr=split(/\s+/,$fp);
    if(exists $listmtsnv{$arr[0]})
    {
      $listmtsnv{$arr[0]}++; 
    }
    else
    {
      $listmtsnv{$arr[0]}=1; 
    }
    undef(@arr);
}
close(FP); 

open(FP,"../datasets/Cancer_mtGenome_analysis_NatGenet2020/TCMA-MutationINDEL.tsv") or die; 
while($fp=<FP>)
{
    chomp($fp); 
    if($fp=~/^sample_id/) { next; }
    @arr=split(/\s+/,$fp);
    if(exists $listmtindel{$arr[0]})
    {
      $listmtindel{$arr[0]}++; 
    }
    else
    {
      $listmtindel{$arr[0]}=1; 
    }
    undef(@arr);
}
close(FP); 

open(FP,"../datasets/PCAWG/donors_and_biospecimens/pcawg_sample_sheet.tsv") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^donor/) { next; }
   @arr=split(/\t/,$fp);
   if($arr[10]!~/Normal/)
   {
     @we=split(/\s+/,$arr[10]); 
     while($we[0]=~/\-/)
     {
        $we[0]=~s/\-//;
     }
     while($we[0]=~/\s+/)
     {
        $we[0]=~s/\s+/\_/;
     }
     $samp=$we[0]; 
     undef(@we); 

     $samptodon{$arr[5]}=$arr[3];   ## $arr[5] - aliquot id; $arr[3] - Donor ID 
     ##print "$arr[3] $arr[10] $arr[11]\n"; 
     $id = $arr[3]; 

     if(exists $sampinfo{$id}) 
     {
	if($arr[11] eq "WGS") 
	{
	   $sampinfo{$id}[0]=$samp;
	}
	elsif($arr[11] eq "RNA-Seq") 
	{
	   $sampinfo{$id}[1]=$samp;
	}
     }
     else
     {
	$sampinfo{$id}[0]="NA"; 
	$sampinfo{$id}[1]="NA"; 
	
	if($arr[11] eq "WGS") 
	{
	   $sampinfo{$id}[0]=$samp;
	}
	elsif($arr[11] eq "RNA-Seq") 
	{
	   $sampinfo{$id}[1]=$samp;
	}
     }
     undef $samp; 
     undef $id; 
   }
   undef(@arr); 
}
close(FP); 

open(FP,"../datasets/PCAWG/consensus_cnv/consensus.20170217.purity.ploidy.txt") or die; 
while($fp=<FP>)
{
   chomp($fp);
   if($fp=~/^samplename/) { next; }
   @arr=split(/\s+/,$fp);
   if(exists $samptodon{$arr[5]}) 
   {
     $dnid=$samptodon{$arr[0]}; 
   }
   else { $dnid = "XXXX"; }   
   $listcnv{$dnid}[0]=$arr[1]; ## Purity
   $listcnv{$dnid}[1]=$arr[2]; ## Ploidy
   $listcnv{$dnid}[2]=$arr[4]; ## wgd or no-wgd
   undef(@arr); 
}
close(FP); 
undef %samptodon; 

## Clinical data 
#

open(FP,"../datasets/PCAWG/Clinical_and_Histology/tumour_subtype_consolidation_map.txt") or die; 
while($fp=<FP>)
{
   chomp($fp);
   @arr=split(/\t/,$fp);
   if($fp=~/^\#/) 
   {
     #for($i=0;$i<$no;$i++)
     #{
     #	print "$i $arr[$i]\n"; 
     #}
     for($i=11;$i<23;$i++)
     {
        $head[$i-11]=$arr[$i];	
     }
     $tot=$i-11; 
     next; 
   }

   $id=$arr[9]; 
   $actmap{$arr[9]} = $arr[12]."*$arr[8]"; 
   #if($arr[8] eq "DO52621") { print "$arr[9]:$arr[12]:$arr[8]:\n"; }

   for($i=11;$i<23;$i++)
   {
      if($arr[$i] eq "") { $arr[$i]="NA"; }
      while($arr[$i]=~/\s+/) 
      {
	  $arr[$i]=~s/\s+/\_/;
      }
      $listclin{$id}[$i-11]=$arr[$i];
   }
   undef $id; 
   undef(@arr);
}
close(FP);


print "Reading Nuclear genome data ...\n"; 

open(FP,"../datasets/PCAWG_filtered/final_consensus_passonly.snv_mnv_indel.icgc.public_duplicate_filt.maf") or die; 
while($fp=<FP>)
{
  chomp($fp); 
  chop($fp); 
  if($fp=~/^Hugo/) { next; }
  @arr=split(/\s+/,$fp); 
  $no=@arr;
  $canc=$arr[$no-2]; 
  $id=$arr[$no-1];
  #if($arr[$no-1] eq "DO52621") { print "YES :$id:\n"; } 

  if(exists $listnucmut{$id})
  {
    if($arr[6] eq "SNP" || $arr[6] eq "DNP") 
    {
	$listnucmut{$id}[0] += length($arr[7]);  
    }
    if($arr[6] eq "INS" || $arr[6] eq "DEL") 
    {
	if($arr[6] eq "INS") 
	{ 
	   $lt=length($arr[9]); 	
	   $listnucmut{$id}[1]+=$lt;  
	   undef $lt; 
        }
	if($arr[6] eq "DEL") 
	{  
	   $lt=length($arr[7]); 	
	   $listnucmut{$id}[1]+=$lt;
	   undef $lt; 
        }  
    }
  }
  else
  { 
    $listnucmut{$id}[0]=0;
    $listnucmut{$id}[1]=0;

    if($arr[6] eq "SNP" || $arr[6] eq "DNP") 
    {
	$listnucmut{$id}[0] += length($arr[7]);
    }
    if($arr[6] eq "INS" || $arr[6] eq "DEL") 
    {
	if($arr[6] eq "INS") 
	{ 
	   $lt=length($arr[9]); 	
	   $listnucmut{$id}[1] += $lt;  
	   undef $lt;
        }
	if($arr[6] eq "DEL") 
	{  
	   $lt=length($arr[7]); 	
	   $listnucmut{$id}[1] += $lt;
	   undef $lt;
        }  
    }
  }
  undef $canc; 
  undef $id; 
  undef $no; 
  undef(@arr);
}
close(FP); 


print "Writing results ...\n";

open(FP,"../datasets/PCAWG/Clinical_and_Histology/pcawg_donor_clinical_August2016_v9.txt") or die; 
open(WR,">../results_June2025/PCAWG_MitoData_collected.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\t/,$fp); 
   if($fp=~/^\#/) 
   {
	print WR "TumorType_WGS\tTumorType_RNASeq\t"; 
	for($i=1;$i<12;$i++)
	{
	  if($arr[$i] eq "project_code") 
	  {
	     $arr[$i] = "CancerType" 
	  }
	  if($arr[$i] eq "donor_age_at_diagnosis")
	  {
		$arr[$i].="_years";
	  }
	  if($arr[$i] eq "donor_survival_time")
	  {
		$arr[$i].="_days";
	  }
	  print WR "$arr[$i]\t";
        }
	undef(@arr);
	for($k=0;$k<$tot;$k++)
	{
	   print WR "$head[$k]\t"; 
	}
	print WR "MT_CopyNumber\tMT_SNV\tMT_INDEL\tTumor_Purity\tTumor_Ploidy\tTumor_wgd_status\tNUC_SNV\tNUC_INDEL\tNormal_Tissue\tNormal_Blood\tNormal_Other\tNormal_AvgNuclearCov\tNormal_MtGenomeCov\tNormal_Mt.Nuclear_CovRatio\tTumorType\tRatio_MTCopyNumber_CancerToNormal\tTotal_Num_Mutations\n";
	next; 
   }

   $tumtyp="NA";
   $patid="NA"; 

   if(exists $sampinfo{$arr[2]})
   {
      print WR "$sampinfo{$arr[2]}[0]\t$sampinfo{$arr[2]}[1]\t";
   }
   else
   {
      print WR "NA\tNA\t"; 
   }

   $id = 'NA'; 
   $canctype = "NA"; 
   if(exists $actmap{$arr[3]}) 
   {
     @sp = split(/\*/,$actmap{$arr[3]});
     $canctype = $sp[0]; 
     if($sp[0] eq "") { $canctype = 'NA'; } 
     $id = $sp[1]; 
     if($sp[1] eq "") { $id= 'NA'; } 
     undef @sp; 
   }


   for($i=1;$i<12;$i++)
   {
      if($arr[$i] eq "") { $arr[$i]="NA"; }
      while($arr[$i]=~/\s+/) 
      {
	  $arr[$i]=~s/\s+/\_/;
      }

      @wr=split(/\-/,$arr[1]);
      $tumtyp=$wr[0];
      undef(@wr); 
      $patid=$arr[2];

      if($i==1) 
      {
	print WR "$canctype\t";  
      }
      else { print WR "$arr[$i]\t"; } 
   }

   #print "$tumtyp\t$patid\n"; 

   if(exists $listclin{$arr[3]}) 
   {
       for($k=0;$k<$tot;$k++)
       {
	   print WR "$listclin{$arr[3]}[$k]\t"; 
       }
   }
   else
   {
       for($k=0;$k<$tot;$k++)
       {
	   print WR "NA\t"; 
       }
   }

   $mtcn="NA"; 
   $totsum=0; 

   if(exists $listmtcn{$arr[3]}) 
   {
	print WR "$listmtcn{$arr[3]}\t";
	$mtcn=$listmtcn{$arr[3]};
   }
   else { print WR "NA\t"; }
   if(exists $listmtsnv{$arr[3]}) 
   {
	print WR "$listmtsnv{$arr[3]}\t";
	$totsum += $listmtsnv{$arr[3]};
   }
   else { print WR "NA\t"; }
   if(exists $listmtindel{$arr[3]}) 
   {
	print WR "$listmtindel{$arr[3]}\t";
	$totsum += $listmtindel{$arr[3]};
   }
   else { print WR "NA\t"; }
   
   if(exists $listcnv{$arr[2]}) 
   {
	print WR "$listcnv{$arr[2]}[0]\t$listcnv{$arr[2]}[1]\t$listcnv{$arr[2]}[2]\t";
   }
   else { print WR "NA\tNA\tNA\t"; }

   
   if(exists $listnucmut{$id}) 
   {
	print WR "$listnucmut{$id}[0]\t$listnucmut{$id}[1]\t";
	
	$totsum += $listnucmut{$id}[0];
	$totsum += $listnucmut{$id}[1];

   }
   else { print WR "NA\tNA\t"; }

   $nrmcn="NA"; 
   $fl=0;
   if($tumtyp!~/NA/ && $patid!~/NA/)
   {
	#print "$tumtyp $patid\n"; 
	$flname=$tumtyp."_Stats.txt"; 
	if(-e "../results/Stats_NormalSamp_Cov_AlignedReads_WgsWxs/$flname")
	{
	  open(TP,"../results/Stats_NormalSamp_Cov_AlignedReads_WgsWxs/$flname") or die; 
          while($tp=<TP>)
	  {
	     chomp($tp); 
             if($tp=~/^SampleType/) { next; }  
	     @tr=split(/\s+/,$tp); 
	     #print "$tr[1] $patid\n"; 
             if($tr[1] eq $patid)
	     {
	       if($tr[0] eq "Tissue")
	       {
	 	  print WR "Tissue\tNA\tNA\t$tr[2]\t$tr[3]\t$tr[4]\t";
	       }
	       if($tr[0] eq "Blood")
	       {
		  print WR "NA\tBlood\tNA\t$tr[2]\t$tr[3]\t$tr[4]\t";
	       }
	       if($tr[0] eq "Other")
	       {
		  print WR "NA\tNA\tOther\t$tr[2]\t$tr[3]\t$tr[4]\t";
	       } 
	       $nrmcn=$tr[4]; 
               $fl=1; 
	    }
	    undef(@tr); 
	  }
	  close(TP);
	}
	undef $flname; 
   }
   if($fl==0)
   {
      print WR "NA\tNA\tNA\tNA\tNA\tNA\t"; 
   }
   print WR "$tumtyp\t";
   $rat="NA"; 
   if(($mtcn ne "NA") && ($nrmcn ne "NA") && ($nrmcn!=0))
   {
      $rat=sprintf("%0.2f",$mtcn/$nrmcn); 
   }
   print WR "$rat\t";
   print WR "$totsum\n";  
   undef $mtcn; 
   undef $nrmcn;
   undef $rat; 
   undef $tumtyp;
   undef $patid; 
   undef $id;
   undef $canctype; 
   undef(@arr);
}
close(WR);
close(FP);

undef(@head);
undef $tot; 

undef %sampinfo; 
undef %listclin; 
undef %listmtcn;
undef %listmtsnv;
undef %listmtindel;
undef %listcnv;
undef %listnucmut; 
undef %actmap; 
