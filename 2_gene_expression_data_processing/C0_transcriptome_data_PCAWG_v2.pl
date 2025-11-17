#!/usr/bin/perl -w 
#
# Process PCAWG transcriptome data and format 
#

print "Reading PCAWG smaple id ...\n"; 

open(FP,"../datasets/PCAWG/donors_and_biospecimens/pcawg_sample_sheet.tsv") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   if($fp=~/^donor_unique_id/) { next; } 
   @arr=split(/\t/,$fp); 
   $donid=$arr[3];
   $samp=$arr[10]; 
   while($samp=~/\-/) 
   {
      $samp=~s/\-//;
   }
   while($samp=~/\s+/) 
   {
      $samp=~s/\s+/\_/;
   }
   if($arr[11]=~/RNA-Seq/) 
   {
      $alqlist{$arr[5]}=$donid.".$samp"; 
      #print "$donid $samp $arr[11]\n"; 
   }
   undef $donid; 
   undef $samp; 
   undef(@arr); 
}
close(FP);

print "Reading gene-transcript mapping ...\n"; 

open(FP,"../datasets/Homo_sapiens_mapping.txt") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\t/,$fp); 
   $gnname=$arr[3];
   $trname=$arr[1];

   $enslist{$trname}=$gnname; 
   ##print "$arr[2] $gnname $trname\n"; 
   undef $trname; 
   undef $gnname; 
   undef(@arr); 
}
close(FP); 

print "Processing PCAWG RNAseq data ...\n"; 

open(FP,"../datasets/PCAWG/transcriptome/transcript_expression/pcawg.rnaseq.transcript.expr.counts.tsv") or die; 
open(WR,">../datasets/PCAWG/transcriptome/transcript_expression/Processed2_pcawg.rnaseq.transcript.expr.counts.tsv") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\s+/,$fp); 
   $no=@arr;

   if($fp=~/^Feature/)
   {
	print WR "Gene\t"; 
	for($i=1;$i<$no;$i++)
	{
	   if(exists $alqlist{$arr[$i]})
	   {
		print WR "$alqlist{$arr[$i]}\t";
	   }
	   else
	   {
		print WR "$arr[$i]\t"; 
	   }
	}	
	print WR "\n"; 
	undef $no;
	undef(@arr); 
	next; 
   }    
   @qw=split(/\./,$arr[0]);
   $trn=$qw[0]; 
   undef(@qw);

   $ttname = $arr[0]; 

   if(exists $enslist{$trn})
   {
	$ttname = $enslist{$trn};  
   }
   if(exists $exprlist{$ttname}) 
   {
     for($i=1;$i<$no;$i++)
     {
       $exprlist{$ttname}[1][$i-1] += $arr[$i];    
     }
   }
   else
   {
     for($i=1;$i<$no;$i++) 
     {
       $exprlist{$ttname}[1][$i-1] = $arr[$i]; 
     }
     $exprlist{$ttname}[0] = $no-1; 
   }

   undef $trn; 
   undef $no; 
   undef(@arr);    
}
close(FP); 

foreach $key (sort keys %exprlist) 
{
   if($key!~/^ENST/) 
   {
      print WR "$key\t"; 
      for($i=0;$i<$exprlist{$key}[0];$i++) 
      {
	 print WR "$exprlist{$key}[1][$i]\t"; 
      }
      print WR "\n"; 
   }
}
undef %exprlist; 

undef %enslist; 
undef %alqlist; 
close(WR); 
