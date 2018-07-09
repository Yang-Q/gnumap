#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

&blast2sam;

sub blast2sam {
  my %opts = ();
  getopts('s', \%opts);
  die("Usage: blast2sam.pl <in.blastn>\n") if (-t STDIN && @ARGV == 0);
  my ($qlen, $slen, $q, $s, $qbeg, $qend, $sbeg, $send, @sam, @cigar, @cmaux, $show_seq, $x);
  $show_seq = defined($opts{s});
  @sam = (); @sam[0,4,6..8,10] = ('', 255, '*', 0, 0, '*');
	
	while (<>) {
		if (/(\S+)\s+total letters/) {
			$slen = $1; $slen =~ s/,//g;
			last;
		}
	}
	
  while ($line=<>) {
    if (@cigar && (/^Query=/ || /Score =.*bits.*Expect/)) { # find a new entry! so wrapup prev matching...
      &blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend);
      @cigar = ();
    }
    if ($line =~ /^Query=\s(\S+)/) { #1 read
      $sam[0] = $1;
		} elsif ($line =~ /Length\s=\s(\d+)/) { #2 read length
      $qlen = $1;
		} elsif ($line =~ /^>(\S+)/) { #3 target seq
      $sam[2] = $1;
    } elsif ($line =~ /Score\s+=\s+(\S+) bits.+Expect(\(\d+\))?\s+=\s+(\S+)/) { #4. the start of an alignment block
      my ($as, $ev) = (int($1 + .499), $3);
      $ev = "1$ev" if ($ev =~ /^e/);
      @sam[1,3,9,11,12] = (0, 0, '', "AS:i:$as", "EV:Z:$ev");
      @cigar = (); $qbeg = 0;
      @cmaux = (0, 0, 0, '');
    } elsif ($line =~ /Strand=(\S+)\/(\S+)/) { #5. mapped strand
      $sam[1] |= 0x10 if ($2 eq 'Minus'); # set 2nd bit to 1
    } elsif ($line =~ /Query\s+(\d+)\s*(\S+)\s+(\d+)/) { #6. st_bp, querySeq, end_bp
      $q = $2; #querySeq
      unless ($qbeg) {
        $qbeg = $1;
        push(@cigar, ($1-1) . "H") if ($1 > 1);
      }
      $qend = $3;
      if ($show_seq) {
        $x = $q;
        $x =~ s/-//g; $sam[9] .= $x;
      }
    } elsif ($line =~ /Sbjct\:\s(\d+)\s*(\S+)\s(\d+)/) { #7. st_bp, targetSeq, end_bp
			$s = $2; #targetSeq
      if ($sam[1] & 0x10) { #determine mapLoc
        $sam[3] = $3;
				$sbeg=$3;
				$send=$1;
      } else {
        $sam[3] = $1 unless ($sam[3]);
				$sbeg=$1;
				$send=$3;
      }
      &aln2cm($qbeg,$qend,$sbeg,$send,\@cigar, $q, $s, \@cmaux);
    }
  }
  &blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend); #don't forget the last wrapup!!!
}

sub blast_print_sam {
  my ($sam, $cigar, $cmaux, $qrest) = @_;
	
  push(@$cigar, $cmaux->[1] . substr("MDI", $cmaux->[0], 1));
  push(@$cigar, $qrest . 'H') if ($qrest);
  if ($sam->[1] & 0x10) { #do reverse compli. on cigar and consSeq too!
    @$cigar = reverse(@$cigar);
    $sam->[9] = reverse($sam->[9]);
    $sam->[9] =~ tr/atgcrymkswATGCRYMKSW/tacgyrkmswTACGYRKMSW/;
  }
  $sam->[9] = '*' if (!$sam->[9]);
  $sam->[5] = join('', @$cigar);
  print join("\t", @$sam), "\n";
}

sub aln2cm {
  my ($qbeg,$qend,$sbeg,$send, $cigar, $q, $s, $cmaux) = @_;
	
	my ($op,$i);
	my $l = $qend-$qbeg+1;
	if ($l==($send-$sbeg+1)) {
		push(@$cigar,"${l}M");
		return;
	}

  for ($i = 0; $i < $l; ++$i) {
    # set $op
    if (substr($q, $i, 1) eq '-') { $op = 2; }
    elsif (substr($s, $i, 1) eq '-') { $op = 1; }
    else { $op = 0; }
    # for CIGAR
    if ($cmaux->[0] == $op) { #keep growing...
      ++$cmaux->[1];
    } else { #wrapup prev a_transcript
      push(@$cigar, $cmaux->[1] . substr("MDI", $cmaux->[0], 1));
      $cmaux->[0] = $op; #assign new aliTranscript
			$cmaux->[1] = 1;
    }
  }
}
