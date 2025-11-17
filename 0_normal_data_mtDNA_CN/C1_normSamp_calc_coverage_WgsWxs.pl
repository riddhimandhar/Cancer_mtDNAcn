#!/usr/bin/perl -w 

## Genome coverage from aligned reads 
## Store patient-wise with Donor ID 
## Only for Normal/Healthy samples

use LWP::Simple;
use Parallel::ForkManager;

system ("ls /mnt/icgc_alignedReadsWGSWXS_collab/*/*.bam > TMP_alignedReads_WgsWxs.txt"); 

open(MF,"../manifest_files/manifest_alignedReads_WGS_WXS/manifest.collaboratory.1605876429526.tsv") or die; 
while($mf=<MF>)
{
  chomp($mf); 
  if($mf=~/^repo/) { next; }
  @arr=split(/\t/,$mf);
  $donid{$arr[4]}=$arr[9]."_$arr[8]";
  undef(@arr); 
}
close(MF);

open(FP,"../datasets/repository_1615486599.tsv") or die; 
open(WR,">../datasets/ICGC_repos_summary.txt") or die; 
while($fp=<FP>)
{
   chomp($fp);
   @arr=split(/\t/,$fp); 
   if($fp=~/^Access/) 
   { 
      print WR "$arr[3]\t"; 
      for($i=4;$i<8;$i++)
      {
	 if($arr[$i] eq "") { $arr[$i]="NA"; } 
         while($arr[$i]=~/\s+/)
         {
           $arr[$i]=~s/\s+/\_/; 
         }
	 ## print "$i $arr[$i]\n"; 
	 print WR "$arr[$i]\t";
      }
      print WR "$arr[9]\t$arr[13]\t";
      print WR "\n";
      $cnt=$i-4; 
      undef(@arr);
      next; 
   }
   print WR "$arr[3]\t"; 

   for($i=4;$i<8;$i++)
   {
     if($arr[$i] eq "") { $arr[$i]="NA"; }
     while($arr[$i]=~/\-/)
     {
        $arr[$i]=~s/\-//; 
     }
     while($arr[$i]=~/\s+/)
     {
        $arr[$i]=~s/\s+/\_/; 
     }
     $fllist{$arr[3]}[$i-4]=$arr[$i];
     print WR "$arr[$i]\t"; 
   }  
   print WR "$arr[9]\t$arr[13]\t";
   print WR "\n";
   undef(@arr);
}
close(FP); 
undef $cnt; 

$i=0;
open(TM,"TMP_alignedReads_WgsWxs.txt") or die;
while($tm=<TM>)
{
  chomp($tm);
  @arr=split(/\//,$tm);
  $elem=$arr[4];
  undef(@arr);

  if($elem=~/^mini/) 
  { 
     undef $elem; 
     next; 
  } 

  if(exists $fllist{$elem} && ($fllist{$elem}[2]=~/Tumor/ || $fllist{$elem}[2]=~/tumor/ || $fllist{$elem}[2]=~/Tumour/ || $fllist{$elem}[2]=~/tumour/))
  {
     ##print "$elem $fllist{$elem}[2]\n"; 
     undef $elem; 
     next;
  }

  if(exists $donid{$elem}) { $din=$donid{$elem}; }
  else { $din=$elem; }

  $wrfile="../results/NormalSamp_Cov_AlignedReads_WgsWxs/$din";
  @ak=split(/\./,$elem);
  $nn=@ak;
  $exp=$ak[0];
  for($m=1;$m<$nn-1;$m++)
  {
     $exp.=".$ak[$m]";
  }
  undef $nn;
  undef(@ak);
  $din.="_$exp";
  undef $exp;
  $wrfile="../results/NormalSamp_Cov_AlignedReads_WgsWxs/$din";
  
  if(exists $fllist{$elem})
  {
     $wrfinal=$wrfile."_$fllist{$elem}[2]";
  }
  $wrfinal=$wrfinal."_Summary";

  if(-e $wrfinal)
  { 
     system "wc -l  $wrfinal > TLS"; 
     open(TLS,"TLS") or die;
     while($tls=<TLS>)
     {
	chomp($tls); 
        @tt=split(/\s+/,$tls);
        if($tt[0]==0)
        { 
	   ##print "YES\n";
	   system "rm $wrfinal";
           $lkarr[$i]=$tm;
           $i++;
        }
        undef(@tt);
     }
     close(TLS);
     system "rm TLS"; 
  }
  else
  {
    $lkarr[$i]=$tm;
    $i++;
  }

  undef $elem;
  undef $din; 
  undef $wrfile;
  undef $wrfinal; 
}
close(TM); 

print "Done checking list ....\n No. of files to be analyzed : $i\n"; 
## Results folder 

my $pm = Parallel::ForkManager->new(20);

 LINKS:
 foreach $full (@lkarr) 
 {
   $pm->start and next LINKS; # do the fork
   @arr=split(/\//,$full); 
   $elem=$arr[4];
   undef(@arr); 

   if(exists $donid{$elem}) { $din=$donid{$elem}; }
   else { $din=$elem; }

   #$wrfile="../results/Cov_AlignedReads_WgsWxs/$din";
   #if(-e $wrfile) 
   #{ 
     @ak=split(/\./,$elem); 
     $nn=@ak;
     $exp=$ak[0];
     for($i=1;$i<$nn-1;$i++)
     {
       $exp.=".$ak[$i]";
     }
     undef $nn;
     undef(@ak);

     $din.="_$exp";
     undef $exp; 
     $wrfile="../results/NormalSamp_Cov_AlignedReads_WgsWxs/$din"; 
   #}

   print "$elem\t$wrfile\n";
   $locflag=0;
   $lcc=0;

   while($locflag==0 && $lcc<=3)
   {
     $lcc++;
     ##print "  Iteration : $lcc\n";
     system ("bedtools genomecov -bga -ibam $full > $wrfile");

     open(FP,$wrfile) or die;
     while($fp=<FP>)
     {
       chomp($fp);
       @arr=split(/\s+/,$fp);
       if(exists $chr{$arr[0]})
       {
          $chr{$arr[0]}[0]+=($arr[3]*($arr[2]-$arr[1]));
          if($arr[2]>$chr{$arr[0]}[1])
          {
            $chr{$arr[0]}[1]=$arr[2];
          }
       }
       else
       {
          $chr{$arr[0]}[0]=$arr[3]*($arr[2]-$arr[1]);
          $chr{$arr[0]}[1]=$arr[2];
       }
       undef(@arr);
     }
     close(FP);

     ## Insert a check and if required repeat 

     $clc=0;
     for($cj=1;$cj<=22;$cj++)
     {
       if(exists $chr{$cj})
       {
         if($chr{$cj}[0]==0)
         {
           $clc++; 
         }
       }
     }

     if($clc==0) { $locflag=1; }
     ##print "$clc $locflag\n";

     if($locflag==1)
     {
       if(exists $fllist{$elem})
       {
         $wrfinal=$wrfile."_$fllist{$elem}[2]";
       }
       $wrfinal=$wrfinal."_Summary";

       open(RW,">$wrfinal") or die;
       foreach $key (sort keys %chr)
       {
         ##print "$key $chr{$key}[0] $chr{$key}[1]\n";
         $avcov=sprintf("%0.2f",$chr{$key}[0]/($chr{$key}[1]+1));
         print RW "$key\t$avcov\n";
         undef $avcov;
       }
       close(RW);
     }
     
     undef %chr; 
     undef $clc; 
     system "rm $wrfile"; 
   }

   undef $lcc;
   undef $locflag; 
   undef %chr; 
   undef $wrfinal;
   undef $wrfile;
   undef $din; 
   $pm->finish; # do the exit in the child process
 }

$pm->wait_all_children;

undef(@lkarr);

undef %donid; 
undef %fllist; 
## system "rm TMP_alignedReads_WgsWxs.txt"; 
