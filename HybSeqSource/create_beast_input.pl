#!/lusr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Bio::AlignIO;


#input: alignments in nexus format
#output: xml file for beast

#assumptions:
#no missing taxa in any gene alignments
#all unlinked substitution models (gtr, estimated, gamma+invariant, 4 gamma categories
#all unlinked clock models (strict, estimate all)
#all unlinked tree priors (birth-death species tree prior, piecewise linear and constant root)
#estimated root height: choose rand between ???
#maps taxa "name" to species "sp.name"

sub badInput {
    my $message = "Usage: perl $0
  -alignment-list=<file: list of files containing alignments>
  -constant-pop: assume constant population sizes
  -taxon-groups=<string: list of groups if all the taxa are not together (default), taxa names separated by space, groups separated by comma>
  -length=<number of mcmc steps>
  -sample=<sample every n states>
  -species-prefix=<prefix for mapping genes to species>
  -cleanup=<remove extra files when done>
  -o=<output xml file>
  ";
    print STDERR $message;
    die "\n";
}

GetOptions(
    "alignment-list=s"=>\my $alignment_list,
    "constant-pop"=>\my $constant_pop,
    "taxon-groups=s"=>\my $taxon_groups,
    "length=i"=>\my $length,
    "sample=i"=>\my $sample,
    "species-prefix=s"=>\my $species_prefix,
    "cleanup"=>\my $cleanup,
    "o=s"=>\my $out,
    );

badInput() unless(defined($alignment_list));
badInput() unless(defined($length));
badInput() unless(defined($sample));
badInput() unless(defined($out));

my %taxa_order;
$species_prefix = defined($species_prefix) ? $species_prefix : "sp";

my $beast_xml_template_0 = '<?xml version="1.0" standalone="yes"?>

<beast>

	<!-- The list of taxa analyse (can also include dates/ages).                 -->
	<!-- ntax=n                                                                 -->
	<taxa id="taxa">
';
#		<taxon id="O"/>
#		<taxon id="1"/>
#		<taxon id="2"/>
my $beast_xml_template_1 = '	</taxa>

	<!-- The sequence alignment (each sequence refers to a taxon above).         -->
	<!-- ntax=n nchar=m                                                      -->
';
#	<alignment id="alignment1" dataType="nucleotide">

#		<sequence>
#			<taxon idref="O"/>
#			CTCGATTATGACTTATATTGCGTGGACCAGCTTATGACAGCGTCGTCTCTTACCCACCTGAGTACTCAAACTCTAGCGTGAATTTAAGTGGACAGCAACCATGTCACTTGGCTATCCCATCAGCTGGAA
#		</sequence>
#		<sequence>
#			<taxon idref="1"/>
#			CTCGATTATGACTTATATTGAGTGGACCAGCTTATGACAGCGTCGTCTCTTACCCACCTGAGTACTCAAACTCTAGCGTGAATTTAAGTGGACAGCAACCATGTCACTTGGCTATCCCATCAGCTGGAA
#		</sequence>

#	</alignment>

my $beast_xml_template_2 = '
	<!-- The unique patterns from 1 to end                                       -->
';
#	<!-- npatterns=a                                                            -->
#	<patterns id="Rep13gseqs.seq0.patterns" from="1">
#		<alignment idref="alignment1"/>
#	</patterns>

#	<!-- The unique patterns from 1 to end                                       -->
#	<!-- npatterns=b                                                            -->
#	<patterns id="Rep13gseqs.seq6.patterns" from="1">
#		<alignment idref="alignment2"/>
#	</patterns>

my $beast_xml_template_3 = '
	<!-- A prior assumption that the population size has remained constant       -->

	<!-- throughout the time spanned by the genealogy.                           -->
	<constantSize id="constant" units="substitutions">
		<populationSize>
';
#			<parameter id="constant.popSize" value="0.0056" lower="0.0" upper="Infinity"/>
my $beast_xml_template_4 = '
		</populationSize>
	</constantSize>
';

#	<!-- Generate a random starting tree under the coalescent process            -->
#	<coalescentTree id="Rep13gseqs.seq0.startingTree" rootHeight="0.0067">
#		<taxa idref="taxa"/>
#		<constantSize idref="constant"/>
#	</coalescentTree>

#	<!-- Generate a random starting tree under the coalescent process            -->
#	<coalescentTree id="Rep13gseqs.seq6.startingTree" rootHeight="0.0046">
#		<taxa idref="taxa"/>
#		<constantSize idref="constant"/>
#	</coalescentTree>

#	<!-- Generate a tree model                                                   -->
#	<treeModel id="Rep13gseqs.seq0.treeModel">
#		<coalescentTree idref="Rep13gseqs.seq0.startingTree"/>
#		<rootHeight>
#			<parameter id="Rep13gseqs.seq0.treeModel.rootHeight"/>
#		</rootHeight>
#		<nodeHeights internalNodes="true">
#			<parameter id="Rep13gseqs.seq0.treeModel.internalNodeHeights"/>
#		</nodeHeights>
#		<nodeHeights internalNodes="true" rootNode="true">
#			<parameter id="Rep13gseqs.seq0.treeModel.allInternalNodeHeights"/>
#		</nodeHeights>
#	</treeModel>

#	<!-- The strict clock (Uniform rates across branches)                        -->
#	<strictClockBranchRates id="Rep13gseqs.seq0.branchRates">
#		<rate>
#			<parameter id="Rep13gseqs.seq0.clock.rate" value="1.0" lower="0.0" upper="Infinity"/>
#		</rate>
#	</strictClockBranchRates>

#	<!-- The general time reversible (GTR) substitution model                    -->
#	<gtrModel id="Rep13gseqs.seq0.gtr">

my $beast_xml_template_5 = '
		<frequencies>
			<frequencyModel dataType="nucleotide">
				<frequencies>
';
#					<parameter id="Rep13gseqs.seq0.frequencies" value="0.25 0.25 0.25 0.25"/>
my $beast_xml_template_6 = '
				</frequencies>
			</frequencyModel>
		</frequencies>
';

#		<rateAC>
#			<parameter id="Rep13gseqs.seq0.ac" value="1.0" lower="0.0" upper="Infinity"/>
#		</rateAC>
#		<rateAG>
#			<parameter id="Rep13gseqs.seq0.ag" value="1.0" lower="0.0" upper="Infinity"/>
#		</rateAG>
#		<rateAT>
#			<parameter id="Rep13gseqs.seq0.at" value="1.0" lower="0.0" upper="Infinity"/>
#		</rateAT>
#		<rateCG>
#			<parameter id="Rep13gseqs.seq0.cg" value="1.0" lower="0.0" upper="Infinity"/>
#		</rateCG>
#		<rateGT>
#			<parameter id="Rep13gseqs.seq0.gt" value="1.0" lower="0.0" upper="Infinity"/>
#		</rateGT>
#	</gtrModel>

#	<!-- site model                                                              -->
#	<siteModel id="Rep13gseqs.seq0.siteModel">
#		<substitutionModel>
#			<gtrModel idref="Rep13gseqs.seq0.gtr"/>
#		</substitutionModel>
#		<gammaShape gammaCategories="4">
#			<parameter id="Rep13gseqs.seq0.alpha" value="0.5" lower="0.0" upper="1000.0"/>
#		</gammaShape>
#		<proportionInvariant>
#			<parameter id="Rep13gseqs.seq0.pInv" value="0.5" lower="0.0" upper="1.0"/>
#		</proportionInvariant>
#	</siteModel>

#	<treeLikelihood id="Rep13gseqs.seq0.treeLikelihood" useAmbiguities="false">
#		<patterns idref="Rep13gseqs.seq0.patterns"/>
#		<treeModel idref="Rep13gseqs.seq0.treeModel"/>
#		<siteModel idref="Rep13gseqs.seq0.siteModel"/>
#		<strictClockBranchRates idref="Rep13gseqs.seq0.branchRates"/>
#	</treeLikelihood>

#	<!-- Species definition: binds taxa, species and gene trees                  -->

#	<species id="species">
#		<sp id="outgroup">
#			<taxon idref="O"/>
#		</sp>
#		<sp id="s1">
#			<taxon idref="1"/>
#		</sp>
#
#		...
#
#
#		<!-- Collection of Gene Trees                                                -->
#		<geneTrees id="geneTrees">
#			<treeModel idref="Rep13gseqs.seq0.treeModel"/>
#			<treeModel idref="Rep13gseqs.seq6.treeModel"/>
#		</geneTrees>
#	</species>

my $beast_xml_template_7 = '
	<!-- Species Tree: Provides Per branch demographic function                  -->
	<speciesTree id="sptree" constantRoot="true">
		<species idref="species"/>
';
my $beast_xml_template_7_constant_pop = '
	<!-- Species Tree: Provides Per branch demographic function                  -->
	<speciesTree id="sptree" constantPopulation="true">
		<species idref="species"/>
';
#		<sppSplitPopulations value="0.0056">
my $beast_xml_template_8 = '
			<parameter id="speciesTree.splitPopSize"/>
		</sppSplitPopulations>
	</speciesTree>

	<!-- Species Tree: tree prior                                                -->

	<!-- Species Tree: Birth Death Model                                         -->
	<birthDeathModel id="birthDeath" units="substitutions">
		<birthMinusDeathRate>
			<parameter id="species.birthDeath.meanGrowthRate" value="1.0" lower="0.0" upper="Infinity"/>
		</birthMinusDeathRate>
		<relativeDeathRate>
			<parameter id="species.birthDeath.relativeDeathRate" value="0.5" lower="0.0" upper="1.0"/>
		</relativeDeathRate>
	</birthDeathModel>

	<!-- Species Tree: Likelihood of species tree                                -->

	<!-- Species Tree: Birth Death Model                                         -->
	<speciationLikelihood id="speciation.likelihood">
		<model>
			<birthDeathModel idref="birthDeath"/>
		</model>
		<speciesTree>
			<speciesTree idref="sptree"/>
		</speciesTree>
	</speciationLikelihood>

	<!-- Species Tree: tmrcaStatistic                                            -->
	<tmrcaStatistic id="speciesTree.rootHeight" name="speciesTree.rootHeight">
		<speciesTree idref="sptree"/>
		<mrca>
			<taxa>
';
#				<sp idref="outgroup"/>
#				<sp idref="s1"/>
#				<sp idref="s2"/>

my $beast_xml_template_9 = '
			</taxa>
		</mrca>
	</tmrcaStatistic>

	<!-- Species Tree: Coalescent likelihood for gene trees under species tree   -->
	<speciesCoalescent id="species.coalescent">
		<species idref="species"/>
		<speciesTree idref="sptree"/>
	</speciesCoalescent>

	<!-- Species tree prior: gama2 + gamma4                                      -->
	<mixedDistributionLikelihood id="species.popSizesLikelihood">
		<distribution0>
			<gammaDistributionModel>
				<shape>
					2
				</shape>
				<scale>
';
#					<parameter id="species.popMean" value="0.0056" lower="0.0" upper="Infinity"/>
my $beast_xml_template_10 = '
				</scale>
			</gammaDistributionModel>
		</distribution0>
		<distribution1>
			<gammaDistributionModel>
				<shape>
					4
				</shape>
				<scale>
					<parameter idref="species.popMean"/>
				</scale>
			</gammaDistributionModel>
		</distribution1>
		<data>
			<parameter idref="speciesTree.splitPopSize"/>
		</data>
		<indicators>
';
#			<parameter value="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"/>
my $beast_xml_template_10a = '
		</indicators>
	</mixedDistributionLikelihood>

	<!-- Define operators                                                        -->
	<operators id="operators">
';

#		<scaleOperator scaleFactor="0.75" weight="0.1">
#			<parameter idref="Rep13gseqs.seq0.ac"/>
#		</scaleOperator>
#		<scaleOperator scaleFactor="0.75" weight="0.1">
#			<parameter idref="Rep13gseqs.seq0.ag"/>
#		</scaleOperator>
#		<scaleOperator scaleFactor="0.75" weight="0.1">
#			<parameter idref="Rep13gseqs.seq0.at"/>
#		</scaleOperator>
#		<scaleOperator scaleFactor="0.75" weight="0.1">
#			<parameter idref="Rep13gseqs.seq0.cg"/>
#		</scaleOperator>
#		<scaleOperator scaleFactor="0.75" weight="0.1">
#			<parameter idref="Rep13gseqs.seq0.gt"/>
#		</scaleOperator>
#		<deltaExchange delta="0.01" weight="0.1">
#			<parameter idref="Rep13gseqs.seq0.frequencies"/>
#		</deltaExchange>
#		<scaleOperator scaleFactor="0.75" weight="0.1">
#			<parameter idref="Rep13gseqs.seq0.alpha"/>
#		</scaleOperator>
#		<scaleOperator scaleFactor="0.75" weight="0.1">
#			<parameter idref="Rep13gseqs.seq0.pInv"/>
#		</scaleOperator>

#		<scaleOperator scaleFactor="0.75" weight="3">
#			<parameter idref="Rep13gseqs.seq0.clock.rate"/>
#		</scaleOperator>

my $beast_xml_template_11 = '
		<upDownOperator scaleFactor="0.75" weight="30">
			<up>
';
#				<parameter idref="Rep13gseqs.seq0.clock.rate"/>
#				<parameter idref="Rep13gseqs.seq6.clock.rate"/>
my $beast_xml_template_12 = '
				<parameter idref="species.birthDeath.meanGrowthRate"/>
			</up>
			<down>
				<speciesTree idref="sptree"/>
				<parameter idref="species.popMean"/>
				<parameter idref="speciesTree.splitPopSize"/>
';
#				<parameter idref="Rep13gseqs.seq0.treeModel.allInternalNodeHeights"/>
#				<parameter idref="Rep13gseqs.seq6.treeModel.allInternalNodeHeights"/>
my $beast_xml_template_13 = '
			</down>
		</upDownOperator>
';
#		<subtreeSlide size="6.7E-4" gaussian="true" weight="15">
#			<treeModel idref="Rep13gseqs.seq0.treeModel"/>
#		</subtreeSlide>
#		<narrowExchange weight="15">
#			<treeModel idref="Rep13gseqs.seq0.treeModel"/>
#		</narrowExchange>
#		<wideExchange weight="3">
#			<treeModel idref="Rep13gseqs.seq0.treeModel"/>
#		</wideExchange>
#		<wilsonBalding weight="3">
#			<treeModel idref="Rep13gseqs.seq0.treeModel"/>
#		</wilsonBalding>
#		<scaleOperator scaleFactor="0.75" weight="3">
#			<parameter idref="Rep13gseqs.seq0.treeModel.rootHeight"/>
#		</scaleOperator>
#		<uniformOperator weight="30">
#			<parameter idref="Rep13gseqs.seq0.treeModel.internalNodeHeights"/>
#		</uniformOperator>

#		<upDownOperator scaleFactor="0.75" weight="3">
#			<up>
#				<parameter idref="Rep13gseqs.seq0.clock.rate"/>
#			</up>
#			<down>
#				<parameter idref="Rep13gseqs.seq0.treeModel.allInternalNodeHeights"/>
#			</down>
#		</upDownOperator>

my $beast_xml_template_14 = '
		<scaleOperator scaleFactor="0.9" weight="5">
			<parameter idref="species.popMean"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="species.birthDeath.meanGrowthRate"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="species.birthDeath.relativeDeathRate"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.5" weight="94">
			<parameter idref="speciesTree.splitPopSize"/>
		</scaleOperator>
		<nodeReHeight weight="94">
			<species idref="species"/>
			<speciesTree idref="sptree"/>
		</nodeReHeight>
	</operators>

	<!-- Define MCMC                                                             -->
';
#	<mcmc id="mcmc" chainLength="6000000" autoOptimize="true" operatorAnalysis="Rep17.param.ops">
my $beast_xml_template_15 = '
		<posterior id="posterior">
			<prior id="prior">
				<speciesCoalescent idref="species.coalescent"/>
				<mixedDistributionLikelihood idref="species.popSizesLikelihood"/>
				<speciationLikelihood idref="speciation.likelihood"/>
';
#				<gammaPrior shape="0.05" scale="10.0" offset="0.0">
#					<parameter idref="Rep13gseqs.seq0.ac"/>
#				</gammaPrior>
#				<gammaPrior shape="0.05" scale="20.0" offset="0.0">
#					<parameter idref="Rep13gseqs.seq0.ag"/>
#				</gammaPrior>
#				<gammaPrior shape="0.05" scale="10.0" offset="0.0">
#					<parameter idref="Rep13gseqs.seq0.at"/>
#				</gammaPrior>
#				<gammaPrior shape="0.05" scale="10.0" offset="0.0">
#					<parameter idref="Rep13gseqs.seq0.cg"/>
#				</gammaPrior>
#				<gammaPrior shape="0.05" scale="10.0" offset="0.0">
#					<parameter idref="Rep13gseqs.seq0.gt"/>
#				</gammaPrior>

#				<gammaPrior shape="0.1" scale="10.0" offset="0.0">
#					<parameter idref="Rep13gseqs.seq0.clock.rate"/>
#				</gammaPrior>

my $beast_xml_template_16 = '
				<oneOnXPrior>
					<parameter idref="species.popMean"/>
				</oneOnXPrior>
				<oneOnXPrior>
					<parameter idref="species.birthDeath.meanGrowthRate"/>
				</oneOnXPrior>
			</prior>
			<likelihood id="likelihood">
';
#				<treeLikelihood idref="Rep13gseqs.seq0.treeLikelihood"/>
#				<treeLikelihood idref="Rep13gseqs.seq6.treeLikelihood"/>
my $beast_xml_template_17 = '
			</likelihood>
		</posterior>
		<operators idref="operators"/>

		<!-- write log to screen                                                     -->
		<log id="screenLog" logEvery="1000">
			<column label="Posterior" dp="4" width="12">
				<posterior idref="posterior"/>
			</column>
			<column label="Prior" dp="4" width="12">
				<prior idref="prior"/>
			</column>
			<column label="Likelihood" dp="4" width="12">
				<likelihood idref="likelihood"/>
			</column>
			<column label="PopMean" dp="4" width="12">
				<parameter idref="species.popMean"/>
			</column>
';
#			<column label="Rep13gseqs.seq0.rootHeight" sf="6" width="12">
#				<parameter idref="Rep13gseqs.seq0.treeModel.rootHeight"/>
#			</column>
#			<column label="Rep13gseqs.seq6.rootHeight" sf="6" width="12">
#				<parameter idref="Rep13gseqs.seq6.treeModel.rootHeight"/>
#			</column>
#			<column label="Rep13gseqs.seq0.clock.rate" sf="6" width="12">
#				<parameter idref="Rep13gseqs.seq0.clock.rate"/>
#			</column>
#			<column label="Rep13gseqs.seq6.clock.rate" sf="6" width="12">
#				<parameter idref="Rep13gseqs.seq6.clock.rate"/>
#			</column>
my $beast_xml_template_18 = '
		</log>

		<!-- write log to file                                                       -->
';
#		<log id="fileLog" logEvery="600" fileName="Rep17.param.log" overwrite="false">
my $beast_xml_template_19 = '
			<posterior idref="posterior"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<speciesCoalescent idref="species.coalescent"/>
			<mixedDistributionLikelihood idref="species.popSizesLikelihood"/>
			<speciationLikelihood idref="speciation.likelihood"/>
			<parameter idref="species.popMean"/>
			<parameter idref="speciesTree.splitPopSize"/>
			<parameter idref="species.birthDeath.meanGrowthRate"/>
			<parameter idref="species.birthDeath.relativeDeathRate"/>
			<tmrcaStatistic idref="speciesTree.rootHeight"/>
';

#			<parameter idref="Rep13gseqs.seq0.treeModel.rootHeight"/>
#			<parameter idref="Rep13gseqs.seq0.ac"/>
#			<parameter idref="Rep13gseqs.seq0.ag"/>
#			<parameter idref="Rep13gseqs.seq0.at"/>
#			<parameter idref="Rep13gseqs.seq0.cg"/>
#			<parameter idref="Rep13gseqs.seq0.gt"/>
#			<parameter idref="Rep13gseqs.seq0.frequencies"/>
#			<parameter idref="Rep13gseqs.seq0.alpha"/>
#			<parameter idref="Rep13gseqs.seq0.pInv"/>
#			<parameter idref="Rep13gseqs.seq0.clock.rate"/>
#			<treeLikelihood idref="Rep13gseqs.seq0.treeLikelihood"/>
my $beast_xml_template_20 = '
		</log>

		<!-- write tree log to file                                                  -->
';
#		<logTree id="species.treeFileLog" logEvery="600" nexusFormat="true" fileName="Rep17.param.species.trees" sortTranslationTable="true">
my $beast_xml_template_21 = '
			<speciesTree idref="sptree"/>
			<posterior idref="posterior"/>
		</logTree>
';
#		<logTree id="Rep13gseqs.seq0.treeFileLog" logEvery="600" nexusFormat="true" fileName="Rep17.param.Rep13gseqs.seq0.(time).trees" sortTranslationTable="true">
#			<treeModel idref="Rep13gseqs.seq0.treeModel"/>
#			<strictClockBranchRates idref="Rep13gseqs.seq0.branchRates"/>
#			<posterior idref="posterior"/>
#		</logTree>
#		<logTree id="Rep13gseqs.seq0.substTreeFileLog" logEvery="600" nexusFormat="true" fileName="Rep17.param.Rep13gseqs.seq0.(subst).trees" branchLengths="substitutions">
#			<treeModel idref="Rep13gseqs.seq0.treeModel"/>
#			<strictClockBranchRates idref="Rep13gseqs.seq0.branchRates"/>
#		</logTree>
my $beast_xml_template_22 = '
	</mcmc>
	<report>
		<property name="timer">
			<mcmc idref="mcmc"/>
		</property>
	</report>
</beast>
';

##################
### start main ###
my @alignment_list = `cat $alignment_list`;
my @alignment_data;
for(my $i=0; $i<scalar(@alignment_list); $i++) {
    $alignment_data[$i] = {};
    #my ($name, undef, undef) = fileparse($alignment_list[$i], qr/\.[^.]*/);
    my ($name, undef, undef) = fileparse($alignment_list[$i]);
    chomp $name;
    $alignment_data[$i]{'name'} = $name;
    $alignment_data[$i]{'seqs'} = get_sequences($alignment_list[$i]);
}

my $pop = (int(rand(21))+40)/10000; # 0.0040 - 0.0060
my @rootHeight;
for(my $i=0; $i<scalar(@alignment_list); $i++) {
    my $height = (int(rand(28))+46)/10000; # 0.0046 - 0.0073
    $rootHeight[$i] = $height;
}
#to set population and root height for larger taxa, using linear interpolation of two (taxa, pop) pairs that BEAUTi produced: (17, 0.0054) and (100, 0.2)
#this gives m = 0.002345 and b = -0.03446
my $numtaxa = scalar(keys(%{$alignment_data[0]{'seqs'}}));
if($numtaxa > 17) {
    $pop = 0.002345*$numtaxa - 0.03446;
#print STDERR "using pop of $pop\n";
    for(my $i=0; $i<scalar(@alignment_list); $i++) {
	$rootHeight[$i] = $pop;
    }
}

my $beast_xml = $beast_xml_template_0;
foreach my $taxa (sort sort_taxa keys %{$alignment_data[0]{'seqs'}}) {
    $beast_xml .= "                <taxon id=\"$taxa\"/>\n";
}

if(defined($taxon_groups)) {
    my @groups = split(/,/, $taxon_groups);
    my $i=0;
    foreach my $group (@groups) {
	$beast_xml .= "        </taxa>\n";
	$beast_xml .= "        <taxa id=\"group$i\">\n";
	$i++;
	
	my @taxon_list = split(/ /, $group);
	foreach my $taxon (@taxon_list) {
	    if($taxon ne "") {
		$beast_xml .= "                <taxon idref=\"$taxon\"/>\n";
	    }
	}
    }
}

$beast_xml .= $beast_xml_template_1;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $num = $i+1;
    $beast_xml .= "	<alignment id=\"alignment$num\" dataType=\"nucleotide\">\n";
    foreach my $taxa (sort sort_taxa keys %{$alignment_data[$i]{'seqs'}}) {
	$beast_xml .= "		<sequence>\n";
	$beast_xml .= "			<taxon idref=\"$taxa\"/>\n";
	$beast_xml .= "			$alignment_data[$i]{'seqs'}{$taxa}\n";
	$beast_xml .= "		</sequence>\n";
    }
    $beast_xml .= "	</alignment>\n";
}

$beast_xml .= $beast_xml_template_2;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $num = $i+1;
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "	<patterns id=\"$name.patterns\" from=\"1\">\n";
    $beast_xml .= "		<alignment idref=\"alignment$num\"/>\n";
    $beast_xml .= "	</patterns>\n";

}

$beast_xml .= $beast_xml_template_3;
$beast_xml .= "			<parameter id=\"constant.popSize\" value=\"$pop\" lower=\"0.0\" upper=\"Infinity\"/>\n";

$beast_xml .= $beast_xml_template_4;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "	<!-- Generate a random starting tree under the coalescent process            -->\n";
    $beast_xml .= "	<coalescentTree id=\"$name.startingTree\" rootHeight=\"$rootHeight[$i]\">\n";
    $beast_xml .= "		<taxa idref=\"taxa\"/>\n";
    $beast_xml .= "		<constantSize idref=\"constant\"/>\n";
    $beast_xml .= "	</coalescentTree>\n";
}

for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "	<!-- Generate a tree model                                                   -->\n";
    $beast_xml .= "	<treeModel id=\"$name.treeModel\">\n";
    $beast_xml .= "		<coalescentTree idref=\"$name.startingTree\"/>\n";
    $beast_xml .= "		<rootHeight>\n";
    $beast_xml .= "			<parameter id=\"$name.treeModel.rootHeight\"/>\n";
    $beast_xml .= "		</rootHeight>\n";
    $beast_xml .= "		<nodeHeights internalNodes=\"true\">\n";
    $beast_xml .= "			<parameter id=\"$name.treeModel.internalNodeHeights\"/>\n";
    $beast_xml .= "		</nodeHeights>\n";
    $beast_xml .= "		<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n";
    $beast_xml .= "			<parameter id=\"$name.treeModel.allInternalNodeHeights\"/>\n";
    $beast_xml .= "		</nodeHeights>\n";
    $beast_xml .= "	</treeModel>\n";
}

for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "	<!-- The strict clock (Uniform rates across branches)                        -->\n";
    $beast_xml .= "	<strictClockBranchRates id=\"$name.branchRates\">\n";
    $beast_xml .= "		<rate>\n";
    $beast_xml .= "			<parameter id=\"$name.clock.rate\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>\n";
    $beast_xml .= "		</rate>\n";
    $beast_xml .= "	</strictClockBranchRates>\n";
}

for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "	<!-- The general time reversible (GTR) substitution model                    -->\n";
    $beast_xml .= "	<gtrModel id=\"$name.gtr\">\n";
    $beast_xml .= $beast_xml_template_5;
    #if we use estimated or all equal, then these parameters appear
    #if we use empirical, then they change to something else
    #under all equal, these are not used as parameters in MCMC (some settings later on are removed)
    $beast_xml .= "					<parameter id=\"$name.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>\n";
    $beast_xml .= $beast_xml_template_6;
    $beast_xml .= "		<rateAC>\n";
    $beast_xml .= "			<parameter id=\"$name.ac\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>\n";
    $beast_xml .= "		</rateAC>\n";
    $beast_xml .= "		<rateAG>\n";
    $beast_xml .= "			<parameter id=\"$name.ag\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>\n";
    $beast_xml .= "		</rateAG>\n";
    $beast_xml .= "		<rateAT>\n";
    $beast_xml .= "			<parameter id=\"$name.at\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>\n";
    $beast_xml .= "		</rateAT>\n";
    $beast_xml .= "		<rateCG>\n";
    $beast_xml .= "			<parameter id=\"$name.cg\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>\n";
    $beast_xml .= "		</rateCG>\n";
    $beast_xml .= "		<rateGT>\n";
    $beast_xml .= "			<parameter id=\"$name.gt\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>\n";
    $beast_xml .= "		</rateGT>\n";
    $beast_xml .= "	</gtrModel>\n";
    $beast_xml .= "\n";
    $beast_xml .= "	<!-- site model                                                              -->\n";
    $beast_xml .= "	<siteModel id=\"$name.siteModel\">\n";
    $beast_xml .= "		<substitutionModel>\n";
    $beast_xml .= "			<gtrModel idref=\"$name.gtr\"/>\n";
    $beast_xml .= "		</substitutionModel>\n";
    $beast_xml .= "		<gammaShape gammaCategories=\"4\">\n";
    $beast_xml .= "			<parameter id=\"$name.alpha\" value=\"0.5\" lower=\"0.0\" upper=\"1000.0\"/>\n";
    $beast_xml .= "		</gammaShape>\n";
    $beast_xml .= "		<proportionInvariant>\n";
    $beast_xml .= "			<parameter id=\"$name.pInv\" value=\"0.5\" lower=\"0.0\" upper=\"1.0\"/>\n";
    $beast_xml .= "		</proportionInvariant>\n";
    $beast_xml .= "	</siteModel>\n";
}

for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "	<treeLikelihood id=\"$name.treeLikelihood\" useAmbiguities=\"false\">\n";
    $beast_xml .= "		<patterns idref=\"$name.patterns\"/>\n";
    $beast_xml .= "		<treeModel idref=\"$name.treeModel\"/>\n";
    $beast_xml .= "		<siteModel idref=\"$name.siteModel\"/>\n";
    $beast_xml .= "		<strictClockBranchRates idref=\"$name.branchRates\"/>\n";
    $beast_xml .= "	</treeLikelihood>\n";
}

$beast_xml .= "	<!-- Species definition: binds taxa, species and gene trees                  -->\n";
$beast_xml .= "	<species id=\"species\">\n";

foreach my $taxa (sort sort_taxa keys %{$alignment_data[0]{'seqs'}}) {
    $beast_xml .= "		<sp id=\"$species_prefix$taxa\">\n"; #mapping done here
    $beast_xml .= "			<taxon idref=\"$taxa\"/>\n";
    $beast_xml .= "		</sp>\n";
}
$beast_xml .= "		<!-- Collection of Gene Trees                                                -->\n";
$beast_xml .= "		<geneTrees id=\"geneTrees\">\n";
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "			<treeModel idref=\"$name.treeModel\"/>\n";
}
$beast_xml .= "		</geneTrees>\n";
$beast_xml .= "	</species>\n";

if($constant_pop) {
    $beast_xml .= $beast_xml_template_7_constant_pop;
}
else {
    $beast_xml .= $beast_xml_template_7;
}
$beast_xml .= "		<sppSplitPopulations value=\"$pop\">\n";

$beast_xml .= $beast_xml_template_8;
foreach my $taxa (sort sort_taxa keys %{$alignment_data[0]{'seqs'}}) {
    $beast_xml .= "				<sp idref=\"$species_prefix$taxa\"/>\n"; #mapping done here

}

$beast_xml .= $beast_xml_template_9;
$beast_xml .= "					<parameter id=\"species.popMean\" value=\"$pop\" lower=\"0.0\" upper=\"Infinity\"/>\n";

#what is the indicator parameter in here?
#constant population: 2n-1 0's
#piece-wise linear and constant root: n 1's and 2n-2 0's
$beast_xml .= $beast_xml_template_10;
my $pop_param = "			<parameter value=\"";
if($constant_pop) {
    for(my $i=0; $i<2*scalar(keys(%{$alignment_data[0]{'seqs'}})) - 2; $i++) {
	$pop_param .= "0 ";
    }
    $pop_param .= "0\"/>\n"
}
else {
    for(my $i=0; $i<scalar(keys(%{$alignment_data[0]{'seqs'}})); $i++) {
	$pop_param .= "1 ";
    }
    for(my $i=0; $i<2*scalar(keys(%{$alignment_data[0]{'seqs'}})) - 3; $i++) {
	$pop_param .= "0 ";
    }
    $pop_param .= "0\"/>\n"
}
$beast_xml .= $pop_param;

$beast_xml .= $beast_xml_template_10a;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n";
    $beast_xml .= "			<parameter idref=\"$name.ac\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n";
    $beast_xml .= "			<parameter idref=\"$name.ag\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n";
    $beast_xml .= "			<parameter idref=\"$name.at\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n";
    $beast_xml .= "			<parameter idref=\"$name.cg\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n";
    $beast_xml .= "			<parameter idref=\"$name.gt\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
    unless($constant_pop) {
	$beast_xml .= "		<deltaExchange delta=\"0.01\" weight=\"0.1\">\n";
	$beast_xml .= "			<parameter idref=\"$name.frequencies\"/>\n";
	$beast_xml .= "		</deltaExchange>\n";
    }
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n";
    $beast_xml .= "			<parameter idref=\"$name.alpha\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n";
    $beast_xml .= "			<parameter idref=\"$name.pInv\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
}
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"3\">\n";
    $beast_xml .= "			<parameter idref=\"$name.clock.rate\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
}

$beast_xml .= $beast_xml_template_11;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "				<parameter idref=\"$name.clock.rate\"/>\n";
}

$beast_xml .= $beast_xml_template_12;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "				<parameter idref=\"$name.treeModel.allInternalNodeHeights\"/>\n";
}

$beast_xml .= $beast_xml_template_13;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    my $height = ($rootHeight[$i]*1000) . "E-4";
    $beast_xml .= "		<subtreeSlide size=\"$height\" gaussian=\"true\" weight=\"15\">\n";
    $beast_xml .= "			<treeModel idref=\"$name.treeModel\"/>\n";
    $beast_xml .= "		</subtreeSlide>\n";
    $beast_xml .= "		<narrowExchange weight=\"15\">\n";
    $beast_xml .= "			<treeModel idref=\"$name.treeModel\"/>\n";
    $beast_xml .= "		</narrowExchange>\n";
    $beast_xml .= "		<wideExchange weight=\"3\">\n";
    $beast_xml .= "			<treeModel idref=\"$name.treeModel\"/>\n";
    $beast_xml .= "		</wideExchange>\n";
    $beast_xml .= "		<wilsonBalding weight=\"3\">\n";
    $beast_xml .= "			<treeModel idref=\"$name.treeModel\"/>\n";
    $beast_xml .= "		</wilsonBalding>\n";
    $beast_xml .= "		<scaleOperator scaleFactor=\"0.75\" weight=\"3\">\n";
    $beast_xml .= "			<parameter idref=\"$name.treeModel.rootHeight\"/>\n";
    $beast_xml .= "		</scaleOperator>\n";
    $beast_xml .= "		<uniformOperator weight=\"30\">\n";
    $beast_xml .= "			<parameter idref=\"$name.treeModel.internalNodeHeights\"/>\n";
    $beast_xml .= "		</uniformOperator>\n";
}
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "		<upDownOperator scaleFactor=\"0.75\" weight=\"3\">\n";
    $beast_xml .= "			<up>\n";
    $beast_xml .= "				<parameter idref=\"$name.clock.rate\"/>\n";
    $beast_xml .= "			</up>\n";
    $beast_xml .= "			<down>\n";
    $beast_xml .= "				<parameter idref=\"$name.treeModel.allInternalNodeHeights\"/>\n";
    $beast_xml .= "			</down>\n";
    $beast_xml .= "		</upDownOperator>\n";
}

$beast_xml .= $beast_xml_template_14;
$beast_xml .= "	<mcmc id=\"mcmc\" chainLength=\"$length\" autoOptimize=\"true\" operatorAnalysis=\"$out.ops\">\n";
$beast_xml .= $beast_xml_template_15;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "				<gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">\n";
    $beast_xml .= "					<parameter idref=\"$name.ac\"/>\n";
    $beast_xml .= "				</gammaPrior>\n";
    $beast_xml .= "				<gammaPrior shape=\"0.05\" scale=\"20.0\" offset=\"0.0\">\n";
    $beast_xml .= "					<parameter idref=\"$name.ag\"/>\n";
    $beast_xml .= "				</gammaPrior>\n";
    $beast_xml .= "				<gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">\n";
    $beast_xml .= "					<parameter idref=\"$name.at\"/>\n";
    $beast_xml .= "				</gammaPrior>\n";
    $beast_xml .= "				<gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">\n";
    $beast_xml .= "					<parameter idref=\"$name.cg\"/>\n";
    $beast_xml .= "				</gammaPrior>\n";
    $beast_xml .= "				<gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">\n";
    $beast_xml .= "					<parameter idref=\"$name.gt\"/>\n";
    $beast_xml .= "				</gammaPrior>\n";
}
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "				<gammaPrior shape=\"0.1\" scale=\"10.0\" offset=\"0.0\">\n";
    $beast_xml .= "					<parameter idref=\"$name.clock.rate\"/>\n";
    $beast_xml .= "				</gammaPrior>\n";
}

$beast_xml .= $beast_xml_template_16;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "				<treeLikelihood idref=\"$name.treeLikelihood\"/>\n";

}

$beast_xml .= $beast_xml_template_17;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "			<column label=\"$name.rootHeight\" sf=\"6\" width=\"12\">\n";
    $beast_xml .= "				<parameter idref=\"$name.treeModel.rootHeight\"/>\n";
    $beast_xml .= "			</column>\n";
}
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "			<column label=\"$name.clock.rate\" sf=\"6\" width=\"12\">\n";
    $beast_xml .= "				<parameter idref=\"$name.clock.rate\"/>\n";
    $beast_xml .= "			</column>\n";
}

$beast_xml .= $beast_xml_template_18;
$beast_xml .= "		<log id=\"fileLog\" logEvery=\"$sample\" fileName=\"$out.log\" overwrite=\"false\">\n";

$beast_xml .= $beast_xml_template_19;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "			<parameter idref=\"$name.treeModel.rootHeight\"/>\n";
    $beast_xml .= "			<parameter idref=\"$name.ac\"/>\n";
    $beast_xml .= "			<parameter idref=\"$name.ag\"/>\n";
    $beast_xml .= "			<parameter idref=\"$name.at\"/>\n";
    $beast_xml .= "			<parameter idref=\"$name.cg\"/>\n";
    $beast_xml .= "			<parameter idref=\"$name.gt\"/>\n";
    unless($constant_pop) {
	$beast_xml .= "			<parameter idref=\"$name.frequencies\"/>\n";
    }
    $beast_xml .= "			<parameter idref=\"$name.alpha\"/>\n";
    $beast_xml .= "			<parameter idref=\"$name.pInv\"/>\n";
}
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "			<parameter idref=\"$name.clock.rate\"/>\n";
}
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    $beast_xml .= "			<treeLikelihood idref=\"$name.treeLikelihood\"/>\n";
}

$beast_xml .= $beast_xml_template_20;
$beast_xml .= "		<logTree id=\"species.treeFileLog\" logEvery=\"$sample\" nexusFormat=\"true\" fileName=\"$out.species.trees\" sortTranslationTable=\"true\">\n";
$beast_xml .= $beast_xml_template_21;
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    #$beast_xml .= "		<logTree id=\"$name.treeFileLog\" logEvery=\"$sample\" nexusFormat=\"true\" fileName=\"$out.$name.(time).trees\" sortTranslationTable=\"true\">\n";
    unless($cleanup) {
	$beast_xml .= "		<logTree id=\"$name.treeFileLog\" logEvery=\"$sample\" nexusFormat=\"true\" fileName=\"$out.$name.trees\" sortTranslationTable=\"true\">\n";
	$beast_xml .= "			<treeModel idref=\"$name.treeModel\"/>\n";
	$beast_xml .= "			<strictClockBranchRates idref=\"$name.branchRates\"/>\n";
	$beast_xml .= "			<posterior idref=\"posterior\"/>\n";
	$beast_xml .= "		</logTree>\n";
    }
}

#optional output; we don't use this because of space constraints
for(my $i=0; $i<scalar(@alignment_data); $i++) {
    my $name = $alignment_data[$i]{'name'};
    #$beast_xml .= "		<logTree id=\"$name.substTreeFileLog\" logEvery=\"$sample\" nexusFormat=\"true\" fileName=\"$out.$name.(subst).trees\" sortTranslationTable=\"true\">\n";
    #$beast_xml .= "			<treeModel idref=\"$name.treeModel\"/>\n";
    #$beast_xml .= "			<strictClockBranchRates idref=\"$name.branchRates\"/>\n";
    #$beast_xml .= "		</logTree>\n";
}
$beast_xml .= $beast_xml_template_22;

open(OUT, ">", $out) or die "can't open $out: $!\n";
print OUT $beast_xml;
close(OUT);

print STDERR "output at $out\n";
print STDERR "done.\n";

### end main ###
################

sub get_sequences {
    my ($file) = @_;
    my %sequences;
    my $i = 0;

    my $aln_obj = Bio::AlignIO->new(-file => $file);
    my $aln = $aln_obj->next_aln();
    foreach my $seq ($aln->each_seq()) {
	$taxa_order{$seq->display_id()} = $i++;
	$sequences{$seq->display_id()} = $seq->seq();
    }

    #foreach my $k (keys %taxa_order) { print "$k has $taxa_order{$k}\n"; }
    
    return \%sequences;
}

#alphabet before numbers
sub sort_taxa {
    return $taxa_order{$a} <=> $taxa_order{$b};
}
