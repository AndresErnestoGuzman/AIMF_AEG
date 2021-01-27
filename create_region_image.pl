use PDL;
use strict;

use PDL::Image2D;
use PDL::ImageND;


my ($file,$threshold,$fwhm_engordar_pixeles) = @ARGV;

my $i = rfits($file)->squeeze;
my $h = rfitshdr($file);
$threshold = $threshold //4;

if(not defined($fwhm_engordar_pixeles)){
		if(defined($$h{BMAJ})){
			$fwhm_engordar_pixeles = $$h{BMAJ}/$$h{CDELT2};
		}
		else{
			die "Cannot define default fwhm default para engordar.";
		}
}

my $mad = median(abs($i->clump(-1)-median($i->clump(-1))));
my $std = 1.4826 * $mad;
printf("Using std=%.3f mJy/beam and threshold = %.1f * sigma = %.2f mJy\n",1e3*$std,$threshold,$threshold*1e3*$std);
my $label = sprintf('%.1fsigma_%.2fmJy',$threshold,$threshold*1e3*$std);
$threshold *= $std;

$i = $i->setbadtoval(0);
my $mask = ($i->isgood) & ($i > $threshold) ;
#$mask->wfits('auxiliar.fits'); system "ds9 auxiliar.fits";
$mask = filter_small_regions($mask,1.133*$$h{BMIN}*$$h{BMAJ}/$$h{CDELT2}**2);
#$mask->wfits('auxiliar.fits'); system "ds9 auxiliar.fits";

exit print("Not enough pixels in mask. Finished.\n") unless(sum($mask));

my $kernelsize = int($fwhm_engordar_pixeles*4+1);
my $kernel_engordar = exp(-4*log(2)*rvals($kernelsize,$kernelsize,{squared=>1})/$fwhm_engordar_pixeles**2);

$mask = convolveND($mask,$kernel_engordar) > 1/2;

# Inner quarter circle filter


$mask->sethdr($h);
$_ = $file;s/\.fits//;
my $basename = $_;
my $out = $basename."_$label.fits";
my $ctrfilename = "${basename}_$label";
$mask->wfits($out);$_ = $out;s/\.fits$//;

system("casa --nologfile -c 'importfits(fitsimage=\"$out\",imagename=\"$_\", overwrite=True)'");

# my $ord = "ds9 $out -contour nlevels 1 -contour levels 1 -contour smooth 2 -contour yes -contour save $ctrfilename.ctr wcs ICRS -exit";
# print "$ord\n"; system $ord;

# system "uniq $ctrfilename.ctr > $ctrfilename.aux ; mv $ctrfilename.aux $ctrfilename.ctr";
# #system "rm $out";

# ds9_to_crtf("$ctrfilename.ctr");

#########################################################################

sub filter_small_regions{
	my $m = shift; die unless defined($m);
	$m = $m->setnantobad->setbadtoval(0);
	my $limitNpix = shift // 5;
	my $comps = cc8compt($m>0);
	#$comps->wfits('auxiliar.fits'); system "ds9 auxiliar.fits";
	my $nreg = int($comps->clump(-1)->maximum);
	foreach my $idxr (1..$nreg){
		if(sum($m == $idxr)<$limitNpix){
			#print "-- $idxr\n";
			$m *= ($m != $idxr);
			#$m->wfits('auxiliar.fits'); system "ds9 auxiliar.fits";
		}
	}

	return $m;
}

sub ds9_to_crtf{
	
	# CRTF format is pretty inflexible. Be care modifyng the "print OO"
	# statements which create the CASA regions.

	my ($ds9) = @_;
	$ds9 = 'ds9.ctr' unless(defined $ds9);
	system "uniq $ds9 > $ds9.aux ; mv $ds9.aux $ds9";
	$_ = $ds9;s/\.ctr//;
	my $o = "$_.crtf";

	open(DS9,$ds9);open OO,">$o";
	print OO "#CRTFv0 CASA Region Text Format version 0\n";
	my ($level,$ra,$dec,@coords,$wcs);
	my ($first,$last) = (0.0); my ($ra_f,$dec_f,$ra_u,$dec_u); 
	my $in = 0;
	while(<DS9>){
		chomp;
		if(/^global/){
			chomp($wcs = <DS9>); # next line
			die "WCS in ctr file!" unless(($wcs eq 'icrs') or ($wcs eq 'fk5'));
			$wcs = uc($wcs);
		}
		if(/^level=(\S+)/){$level = $1}
		if(/^\(/){
			$in = 1;
			print OO "poly [";
			$first = 1;
			next;
		}
		if($in){
			if(/([0-9]\S*)\s+(\S+)/){
				($ra,$dec) = ($1,$2);
				($ra_f,$dec_f ) = ($ra,$dec) if $first;
				$first = 0 if $first;
				push(@coords,"[${ra}deg, ${dec}deg]");
			}
			if(/^\)/){
				($ra_u,$dec_u) = ($ra,$dec);
				pop @coords if($ra_u == $ra_f and $dec_u==$dec_f);
				print OO join(', ',@coords); # 
				print OO "] coord=$wcs, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=magenta, font=Helvetica, fontsize=11, fontstyle=normal, usetex=false\n";
				$in = 0;
				@coords = ();
				next;
			}
		}
	}
	close DS9;
	close OO;
}
