#!/usr/bin/perl
$thisfilename = "prmtopandpdb2gamess.pl";
#Dan Sindhikara sindhikara@gmail.com
#last update June 3, 2011
sub formatnumber{
	$number=$_[0];
	if ($number>=0) {
		$number=sprintf(" %6.4f",$number);
	} else { $number=sprintf("%6.4f",$number);}
	return($number);
}


print "#Command: "; foreach (@ARGV){print $_;}print "\n";
print "#input format $thisfilename <pdbfile> <prmtoppath> ...\n";
if($#ARGV<1){print "Insufficient arguments!\n";exit(0);}
print "#Warning!! pdb and parmtop must have identical atom type ordering, if not wrong parameters will be assigned\n";
$pdbfilename=@ARGV[0];
$prmtopname=@ARGV[1];
sub parse {
	# need $count @arrayname
       $count=$_[0]; # 
       $line=<$prmtop>;#read FORMAT line
       @splitline=split(/\%FORMAT\(/,$line);
	@splitline=split(/\D/,@splitline[1]);#split on non-digit
	$ncolumns=int(@splitline[0]);
	$colwidth=int(@splitline[1]);
       for $i (1 .. $count) {#
              $column=($i-1)%$ncolumns;
              if ($column==0) {$line=<$prmtop>;}
              @myarray[$i]=(substr($line, int($column*$colwidth),$colwidth));

		}
	return @myarray;
}

open $prmtop,$prmtopname;
#$line=<$prmtop>;print $line;
while ($line=<$prmtop>) {
	@splitline=split("FLAG ",$line);
	@splitline=split(" ",@splitline[1]);
    if (@splitline[0] eq "POINTERS"){
	@pointer= (parse(19));
	for $i (1 .. 19){@pointer[$i]=int(@pointer[$i]);}
	# assign pointers to variables
	$natoms=@pointer[1];
	$ntypes=@pointer[2];
	$nbonh=@pointer[3]; #number of bonds incl H
	$mbona=@pointer[4]; #number of bonds not incl H
	$ntheth=@pointer[5]; #number of angles incl H
	$mtheta=@pointer[6]; #number of angles not incl H
	$nphih=@pointer[7]; #number of dihedrals incl H
	$mphia=@pointer[8]; #number of dihedrals not incl H
	$numbnd=@pointer[16]; #number of unique bonds
	$numang=@pointer[17]; #number of unique angles
	$nptra=@pointer[18]; #number of unique dihedrals
#print "# found $natoms atoms\n";
    }
	if (@splitline[0] eq "ATOM_TYPE_INDEX"){
		@atomtypeindex = parse($natoms);
		for $i (1 .. $natoms){@atomtypeindex[$i]=int(@atomtypeindex[$i]);}
	}
    if (@splitline[0] eq "ATOM_NAME"){@atomname = parse($natoms);}
    if (@splitline[0] eq "CHARGE"){
		@charge = parse($natoms);
    for $i (1 .. $natoms){@charge[$i]=1*(@charge[$i]);}
	}
#       if (@splitline[0] eq "MASS"){
#	@mass = parse($natoms);
#               for $i (1 .. $natoms){@mass[$i]=1*(@mass[$i]);}
#       }
	if (@splitline[0] eq "AMBER_ATOM_TYPE"){
		@amberatomtype = parse($natoms); #string
	}
#if (@splitline[0] eq "BOND_FORCE_CONSTANT"){
#	@bondfc = parse ($numbnd);
#	for $i (1 .. $numbnd){@bondfc[$i]=1*@bondfc[$i];}
#}
#       if (@splitline[0] eq "BOND_EQUIL_VALUE"){
#               @bondev = parse ($numbnd);
#               for $i (1 .. $numbnd){@bondev[$i]=1*@bondev[$i];}
#       }
#       if (@splitline[0] eq "ANGLE_FORCE_CONSTANT"){
#               @anglefc = parse ($numang);
#               for $i (1 .. $numang){@anglefc[$i]=1*@anglefc[$i];}
 #      }
#       if (@splitline[0] eq "ANGLE_EQUIL_VALUE"){
#               @angleev = parse ($numang);
#               for $i (1 .. $numang){@angleev[$i]=1*@angleev[$i];}
#       }
#       if (@splitline[0] eq "DIHEDRAL_FORCE_CONSTANT"){
#               @dihfc = parse ($nptra);
#               for $i (1 .. $nptra){@dihfc[$i]=1*@dihfc[$i];}
#       }
#       if (@splitline[0] eq "DIHEDRAL_PERIODICITY"){
#               @dihpedcty = parse ($nptra);
#               for $i (1 .. $nptra){@dihpedcty[$i]=1*@dihpedcty[$i];}
#       }
#       if (@splitline[0] eq "DIHEDRAL_PHASE"){
#               @dihphase = parse ($nptra);
#               for $i (1 .. $nptra){@dihphase[$i]=1*@dihphase[$i];}
#       }
#       if (@splitline[0] eq "BONDS_INC_HYDROGEN"){
#               @bic = parse ($nbonh*3); # the 3 comes since there are 3 values for each bond
#               for $i (1 .. $nbonh*3){@bic[$i]=int(@bic[$i]);}
#       }
#       if (@splitline[0] eq "BONDS_WITHOUT_HYDROGEN"){
#               @bnic = parse ($mbona*3);
#               for $i (1 .. $mbona*3){@bnic[$i]=int(@bnic[$i]);}
#       }
#       if (@splitline[0] eq "ANGLES_INC_HYDROGEN"){
#               @aic = parse ($ntheth*4); # 4 values per angle 
#               for $i (1 .. $ntheth*4){@aic[$i]=int(@aic[$i]);}
#       }
#if (@splitline[0] eq "ANGLES_WITHOUT_HYDROGEN"){
#	@anic = parse ($mtheta*4);
#	for $i (1 .. $mtheta*4){@anic[$i]=int(@anic[$i]);}
#}
#       if (@splitline[0] eq "DIHEDRALS_INC_HYDROGEN"){
#               @dic = parse ($nphih*5); # 5 values per dihedral 
#               for $i (1 .. $nphih*5){@dic[$i]=int(@dic[$i]);}
#       }
#        if (@splitline[0] eq "DIHEDRALS_WITHOUT_HYDROGEN"){
#                @dnic = parse ($mphia*5); # 5 values per dihedral 
#                for $i (1 .. $mphia*5){@dnic[$i]=int(@dnic[$i]);}
#        }
}
close($prmtop);

print "# Finding and processing non-bonded parameters\n";
system("echo printlennardjones > rdparmin.delete\necho exit >> rdparmin.delete\n");

for $i (1 .. $ntypes){
#	$test=`rdparm $prmtopname < rdparmin.delete | awk \'NR== $i {tot=tot+\$1} END {print tot}\'`;
#	if($test<1){print "Error parmtop does not match pdb atom names\n";exit(0);}
#print "rdparm $prmtopname < rdparmin.delete | awk \'\$1==\$2 && \$2!=\"\" {print \$3}\' | awk \'NR== $i {print}\'\n";
	$T12=`rdparm $prmtopname < rdparmin.delete | grep -A1000 "Lennard-Jones" | awk \'\$1==\$2 && \$2!=\"\" {print \$3}\' | awk \'NR== $i {print}\'`;
    $T6=`rdparm $prmtopname < rdparmin.delete | grep -A1000 "Lennard-Jones" | awk \'\$1==\$2 && \$2!=\"\" {print \$4}\' | awk \'NR== $i {print}\'`;
	#sigma = pow(A/B,1/6), epsilo=B^2/4A
	if($T12==0){@lj1[$i]=0;@lj2[$i]=0;} else {
		@lj1type[$i]=($T12/$T6)**(1/6);
		@lj2type[$i]=($T6*$T6)/(4*$T12);
	}
#printf("@lj1type[$i], @lj2type[$i]\n");


	#	@lj1[$i]=@lj1[$i]*0.8908987;

#		chomp($lj2[$i]);
 #       @amberatomtypelist[@atomtypeindex[$i]]=@amberatomtype[$i];
}
print "# assigning parameters to appropriate atoms\n";

for $i (1 .. $natoms){
	@lj1[$i]=@lj1type[@atomtypeindex[$i]];
	@lj2[$i]=@lj2type[@atomtypeindex[$i]];
}




open $pdbfile,$pdbfilename;
while ($line = <$pdbfile>){
	if ((substr($line,0,6) eq "HETATM") or (substr($line,0,4) eq "ATOM")){
		@atmnum=substr($line,6,5);
		@x[@atmnum]=substr($line,30,8);
                @y[@atmnum]=substr($line,38,8);
                @z[@atmnum]=substr($line,46,8);
#		 @splitline=split(' ',$line);
#		 @x[@splitline[1]]=@splitline[5];
#                @y[@splitline[1]]=@splitline[6];
#                @z[@splitline[1]]=@splitline[7];
	}
}



print "$natoms\n";
for $i (1 .. $natoms){ 
	$nameletter=substr(@atomname[$i],0,4);
	printf("%s  %s   %s   %s   %s   %s   %s\n",$nameletter,formatnumber($lj1[$i]),formatnumber(@lj2[$i]),formatnumber(@charge[$i]/18.2223),formatnumber(@x[$i]),formatnumber(@y[$i]),formatnumber(@z[$i]));

#	printf("%s  1    %5.4f   %5.4f   %*.4f   %5s   %5s   %5s\n",@atomname[$i],$lj1[$i],@lj2[$i],@charge[$i]/18.2223,@x[$i],@y[$i],@z[$i]);


#  CHRG   : the atom charges.  (Divide by 18.2223 to convert to charge

}
