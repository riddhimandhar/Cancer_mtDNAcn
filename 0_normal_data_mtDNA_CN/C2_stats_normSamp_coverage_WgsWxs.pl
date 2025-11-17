#!/usr/bin/perl -w
#

## Statistics on Normal sample WGS data 
#

system "ls ../results/NormalSamp_Cov_AlignedReads_WgsWxs/ > TMSTAT"; 

open(TM,"TMSTAT") or die; 
while($tm=<TM>)
{
   chomp($tm); 
   @qw=split(/\_/,$tm);
   $tumtyp=$qw[0];
   @we=split(/\-/,$tumtyp); 
   $tumsp=$we[0];
   $patid=$qw[1];
   undef(@we);
   undef(@qw); 
   if($tm=~/blood/) { $flag="Blood"; }
   elsif($tm=~/tissue/) { $flag="Tissue"; }
   else { $flag="Other"; }

   if(exists $dup{$patid."*$flag"})
   {
      undef $flag; 
      undef $tumsp; 
      undef $tumtyp; 
      undef $patid; 
      next; 
   }
   $dup{$patid."*$flag"}=1; 
   
   $mtcov=0;
   $avcov=0; $sn=0; 
   open(FP,"../results/NormalSamp_Cov_AlignedReads_WgsWxs/$tm") or die; 
   while($fp=<FP>)
   {
      chomp($fp);
      @arr=split(/\s+/,$fp); 
      if($arr[0] eq "MT")
      {
	  $mtcov=$arr[1]; 
      }
      elsif($arr[0]=~/^\d+/)
      {
          $avcov+=$arr[1];
	  $sn++; 
	  #print "$arr[0]\n"; 
      }
      undef(@arr);
   }
   close(FP); 
   if($sn!=0) 
   { 
	$avcov=sprintf("%0.2f",$avcov/$sn); 
   } 
   ##print "$tumtyp $tumsp $flag $avcov $mtcov\n"; 

   if($avcov==0) { $rat="NA"; }
   else
   {
      $rat=sprintf("%0.2f",$mtcov/$avcov);
   }

   if(exists $normlist{$tumsp}) 
   {
      if($rat!~/NA/)
      {
	$cnt=$alllist{$tumsp}[0];
	$alllist{$tumsp}[1][$cnt]=$flag;
	$alllist{$tumsp}[2][$cnt]=$patid;
	$alllist{$tumsp}[3][$cnt]=$avcov;
	$alllist{$tumsp}[4][$cnt]=$mtcov;
	$alllist{$tumsp}[5][$cnt]=$rat;
	$alllist{$tumsp}[0]++; 
	undef $cnt; 

        $normlist{$tumsp}[0][1]+=$rat;
        $normlist{$tumsp}[0][0]++;

        if($flag eq "Tissue")
        {
          $normlist{$tumsp}[1][1]+=$rat;
          $normlist{$tumsp}[1][0]++;
        } 
        if($flag eq "Blood")
        {
          $normlist{$tumsp}[2][1]+=$rat; 
	  $normlist{$tumsp}[2][0]++;
        } 
        if($flag eq "Other")
        {
          $normlist{$tumsp}[3][1]+=$rat; 
          $normlist{$tumsp}[3][0]++;
	  @we=split(/\_/,$tm);
          $normlist{$tumsp}[4].="$we[4]".".$we[5],";
	  undef(@we);
        }
      }	
   }
   else
   {

      $normlist{$tumsp}[0][0]=0; 
      $normlist{$tumsp}[0][1]=0; 
      $normlist{$tumsp}[1][0]=0; 
      $normlist{$tumsp}[1][1]=0; 
      $normlist{$tumsp}[2][0]=0; 
      $normlist{$tumsp}[2][1]=0; 
      $normlist{$tumsp}[3][0]=0; 
      $normlist{$tumsp}[3][1]=0; 
      $normlist{$tumsp}[4]=""; 

      if($rat!~/NA/)
      {
	$alllist{$tumsp}[1][0]=$flag;
	$alllist{$tumsp}[2][0]=$patid;
	$alllist{$tumsp}[3][0]=$avcov;
	$alllist{$tumsp}[4][0]=$mtcov;
	$alllist{$tumsp}[5][0]=$rat;
	$alllist{$tumsp}[0]=1; 
 
        $normlist{$tumsp}[0][1]+=$rat;
        $normlist{$tumsp}[0][0]++;

        if($flag eq "Tissue")
        {
          $normlist{$tumsp}[1][1]+=$rat;
          $normlist{$tumsp}[1][0]++;
        } 
        if($flag eq "Blood")
        {
          $normlist{$tumsp}[2][1]+=$rat; 
	  $normlist{$tumsp}[2][0]++;
        } 
        if($flag eq "Other")
        {
          $normlist{$tumsp}[3][1]+=$rat; 
          $normlist{$tumsp}[3][0]++;
	  @we=split(/\_/,$tm);
          $normlist{$tumsp}[4]=$we[4].".$we[5],";
	  undef(@we);
        }
      }	
   }

   undef $flag; 
   undef $patid; 
   undef $tumsp; 
   undef $tumtyp; 

   #exit();
}
close(TM); 
undef %dup; 

foreach $key (sort keys %alllist)
{
   $wrfile="../results/Stats_NormalSamp_Cov_AlignedReads_WgsWxs/$key"."_Stats.txt"; 

   open(WR,">$wrfile") or die; 
   print WR "SampleType\tPatientID\tAvgNuclearCov\tMtGenomeCov\tMt.Nuclear_CovRatio\n"; 
   for($i=0;$i<$alllist{$key}[0];$i++)
   {
      print WR "$alllist{$key}[1][$i]\t$alllist{$key}[2][$i]\t$alllist{$key}[3][$i]\t$alllist{$key}[4][$i]\t$alllist{$key}[5][$i]\n";
   } 
   close(WR); 
}

open(WR,">../results/Stats_NormalSamp_Cov_AlignedReads_WgsWxs/Complete_Stats.txt") or die; 
print WR "CancerType\tOverall_Mt.Nuc_CovRat\tOverall_SampleNo\tTissue_Mt.Nuc_CovRat\tTissue_SampleNo\tBlood_Mt.Nuc_CovRat\tBlood_SampleNo\tOther_Mt.Nuc_CovRat\tOther_SampleNo\n"; 
foreach $key (sort keys %normlist)
{
   if($normlist{$key}[0][0]!=0)
   {   
     $normlist{$key}[0][1]=sprintf("%0.2f",$normlist{$key}[0][1]/$normlist{$key}[0][0]);
   }
   if($normlist{$key}[1][0]!=0)
   {   
      $normlist{$key}[1][1]=sprintf("%0.2f",$normlist{$key}[1][1]/$normlist{$key}[1][0]);
   }
   if($normlist{$key}[2][0]!=0)
   {   
     $normlist{$key}[2][1]=sprintf("%0.2f",$normlist{$key}[2][1]/$normlist{$key}[2][0]);
   }
   if($normlist{$key}[3][0]!=0)
   {   
     $normlist{$key}[3][1]=sprintf("%0.2f",$normlist{$key}[3][1]/$normlist{$key}[3][0]);
   }

   ##print "CancerType $key Normal Sample MT/NUC genome ratio: Overall - $normlist{$key}[0][1]\t$normlist{$key}[0][0]\tTissue - $normlist{$key}[1][1]\t$normlist{$key}[1][0]\tBlood - $normlist{$key}[2][1]\t$normlist{$key}[2][0]\tOther - $normlist{$key}[3][1]\t$normlist{$key}[3][0]\n"; 
   print WR "$key $normlist{$key}[0][1]\t$normlist{$key}[0][0]\t$normlist{$key}[1][1]\t$normlist{$key}[1][0]\t$normlist{$key}[2][1]\t$normlist{$key}[2][0]\t$normlist{$key}[3][1]\t$normlist{$key}[3][0]\t$normlist{$key}[4]\n"; 
}
close(WR); 

undef %normlist; 

system "rm TMSTAT"; 
