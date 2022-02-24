use strict;
use warnings;

use IO::Handle;
use IO::File;

my $fh;

if(@ARGV==0){
        $fh = IO::Handle->new();
        $fh->fopen(fileno(STDIN),"r");
} else {
        $fh = IO::File->new();
        $fh->open($ARGV[0],"r");
}

my $id = $ENV{"ID"} || "All";
my $addup = $ENV{"ADDUP"} || "yes";
my $totalReadsAll=0;
my $writtenReadsAll=0;
my $adapterReadsAll=0;
while(my $line=$fh->getline()){
        if($line=~ m/Total reads processed:(.*)/){
                my $totalReads=$1;
                $totalReads=~ s/\s//sg;
                $totalReads=~ s/,//sg;
                $totalReadsAll = $addup eq "yes" ? $totalReads+$totalReadsAll : $totalReads;
        }
        if($line=~ m/Reads with adapters:(.*)\(/){
                my $adapterReads = $1;
                $adapterReads =~ s/\s//sg;
                $adapterReads =~ s/,//sg;
                $adapterReadsAll = $addup eq "yes" ? $adapterReads + $adapterReadsAll : $adapterReads;
        }
        if($line=~ m/Reads written.*:(.*)\(/){
                my $writtenReads = $1;
                $writtenReads =~ s/\s//sg;
                $writtenReads =~ s/,//sg;
                $writtenReadsAll = $addup eq "yes" ? $writtenReads + $writtenReadsAll : $writtenReads;
        }
}
print $id."\t".$totalReadsAll."\t".$writtenReadsAll."\t".$adapterReadsAll."\n";
