# GeneSeqRules : Open-source software for finding biologically relevant sequential rules from longitudinal human microarrays. 

It gathers a set of R, Python and Java-based scripts for performing sequential association rule mining in gene expression temporal data. The scripts and codes available here will allow the user to perform:

	1) Microarray pre-processing (normalisation, probe annotation, analysis and discretization)

	2) Knowledge-extraction (extraction of sequential association rules (gene networks) from discretized gene expression data)

	3) Biological validation of results (integration of output gene networks along with external biological resources such as functional annotation and gene regulation databases (e.g. GO,KEGG and TRRUST))

	4) Visual representation of results (hierarchical edge bundling visualization for the joint representation of gene networks and all accessed biological information)


These processes have been specially designed for Affymetrix array platforms and are organized in seven pipeline STEPs. The instructions for running each part of the code are detailed bellow. For exampling purposes, we have included the raw expression data from the public dataset GSE103766, which can also be downloaded from the next URL (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103766).





#	STEP - 1 	PRE-PREPROCESSING


These lines are instructions for running the "./GeneSeqRules/required_codes/Part_1_Affy_gene_expression_Preprocessing_and_analysis.r" R script, which is intended for the pre-processing and analysis of raw fluorescence (.cel) files. This script performs:

	- The merging of individual .cel files.

	- The visual inspection of raw data and the generation of quality control plots.

	- The RMA normalization process.

	- The annotation process.

	- The differentially expressed gene analysis.


Required Script:

	- "./GeneSeqRules/required_codes/Part_1_Affy_gene_expression_Preprocessing_and_analysis.r"


Instructions:


	(1) Open an Ubuntu-console or a Windows-Command Prompt.

	(2) Run R environment 3.6.0 or higher.

	(3) Install the required packages listed in the script.

	(4) Run the script.


Important Note: The script is adapted to the current example problem. If you change the localization of input required .cel files or their name, you will need to modify the code. Each code line has been carefully annotated with full explanation. Read in detail each code comment before making any change.





#	STEP - 2 	PRE-PREPROCESSING


These lines are instructions for running the "./GeneSeqRules/required_codes/Part_2_gene_expression_discretization_and_sequence_database_generation.r" R script, which is intended for the discretization of absolute gene expression values and the generation of the sequence database. 


Required Script:

	- "./GeneSeqRules/required_codes/Part_2_gene_expression_discretization_and_sequence_database_generation.r"


Instructions:


	(1) Open an Ubuntu-console or a Windows-Command Prompt.

	(2) Run R environment 3.6.0 or higher.

	(3) Install the required packages listed in the script.

	(4) Run the script.


Important Note: The script is adapted to the current example problem. If you change the localization of input required files or their name, you will need to modify the code. The current script is prepared for discretizing a microarray with three temporal records. In the case of analyzing an array dataset with 4 available time records, you will need to adapt the code. Each code line has been carefully annotated with full explanation. Read in detail each code comment before making any change.





#	STEP - 3 	PRE-PREPROCESSING


These lines are instructions for adapting each generated sequence database to the data format required by the CMRules algorithm. Note: the process should be repeated individually in each of the generated sequence databases.


Required Input files:

	- Each generated sequence database. They are available as "./GeneSeqRules/Output_files/Sequence_databases/Group_1/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_1.csv" for group 1 and as "./GeneSeqRules/Output_files/Sequence_databases/Group_2/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_2.csv" for group 2.


Instructions:

	(1) Open each generated sequence database with a spreadsheet manager such as Microsoft Excel and remove row- and column- names (they correspond with the first column and the first row from the dataset respectively).

	(2) Find and replace with an empty space all NAs from the sequence database. If you are doing this modification with Microsoft Excel, please, select the option to move columns to the left when removing NAs.

	(3) Copy and paste the resulting rows and columns into a text editor such as Gedit or Notepad++. Then, replace all "\t" by a single space.

	(4) Save each output sequence database with the name "GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_X_formated" in its corresponding directory. That is to say:

	"./GeneSeqRules/Output_files/Sequence_databases/Group_1/Formated_Sequence_DB/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_1_formated"

	for the file group 1 AND

	"./GeneSeqRules/Output_files/Sequence_databases/Group_2/Formated_Sequence_DB/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_2_formated"

	for the file group 2





#	STEP - 4 	KNOWLEDGE-EXTRACTION
   

These lines are instructions for running the knowledge-extraction stage of our proposed pipeline. Particularly, this code will apply the CMRules algorithm in each generated sequence database. Note: the process should be repeated individually in each of the generated sequence databases. For this part of the pipeline, we employ a third-person software named SPMF. SPMF is an open-source data mining mining library written in Java, specialized in pattern mining (the discovery of patterns in data) (http://www.philippe-fournier-viger.com/spmf/).


Required Input files:

	- Formated sequence databases. They have been generated in the previous step and should be available as "./GeneSeqRules/Output_files/Sequence_databases/Group_1/Formated_Sequence_DB/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_1_formated" for group 1 and as "./GeneSeqRules/Output_files/Sequence_databases/Group_2/Formated_Sequence_DB/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_2_formated" for group 2.


Instructions:

	(1) Open an Ubuntu-console or a Windows-Command Prompt.

	(2) Run (as superuser) the .jar file "./GeneSeqRules/required_codes/spmf.jar".

	(3) Use the GUI of the SPMF library for running the CMRules algorithm in each adapted database. For ensuring an adequate behaviour of the next code scripts and steps, we recommend to chose the next directory and file names for saving each output rule file:

		- Directory: "./GeneSeqRules/Output_files/Output_rules/"
		- Name: "Output_rules_raw_group_1" or "Output_rules_raw_group_2"

		We also recommend to use the next threshold parameters for running the algorithm without memory overload:

		Minsup = 0.6
		Minconf = 0.5
		Max anteceent size = 2
		Max consequent size = 2


	(4) Once the sequential rules have been created with their sequential confidence and sequential support metrics, you can compute the rest of sequential quality metrics by rule (Lift, Conviction, Certainty factor). For that purpose, instead of the ordinal spmf.jar library, you should use our module "calculate_measures_sequential_rules" (available in our modification of the SPMF library "https://github.com/segura2010/SPMF"). For running this  module you can follow the next steps:

		(4.1) Compile our SPMF modification with a java-based build tool such as "Ant" or "Netbeans".

		(4.2) Use netbeans for changing the code of our script "/src/ca/pfv/spmf/tools/calculate_measures_sequential_rules/Main.java". For that purpose, you should change the next code line available in the script:

		"calculator.calculate("/Users/alberto/Desktop/Universidad/Doctorado/NuevasPropuestas/DBs/FIFA_db.txt", "/Users/alberto/Desktop/Universidad/Doctorado/NuevasPropuestas/DBs/resultados/FIFA/cmdeo.txt", "/Users/alberto/Desktop/Universidad/Doctorado/NuevasPropuestas/DBs/resultados/FIFA/cmdeo_measures.txt");"

		Each path in this code line should be respectively replaced with the paths to the next files or directories:

			- The employed input sequence database: "./GeneSeqRules/Output_files/Sequence_databases/Group_1/Formated_Sequence_DB/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_1_formated" for group 1 and "./GeneSeqRules/Output_files/Sequence_databases/Group_2/Formated_Sequence_DB/GENE_EXPRESSION_Sequence_DB_meanFC_GROUP_2_formated" for group 2

			- The obtained output rule files: E.g. "./GeneSeqRules/Output_files/Output_rules/Output_rules_raw_group_1" in the case of group 1

			- The output directory for saving the rules annotated with the rest of quality metrics. We recommend to use the next directory for saving results:

			"/home/augusto/Descargas/IDEA/GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/" with the name "Output_rules_raw_group_1_measures.txt" for group 1 and "Output_rules_raw_group_2_measures.txt"

		(4.3) Re-build the java project for saving code changes and run the script with Netbeans.

Note: If the user want to implement the contrast between two databases when extracting sequential association rules (as we indicate in the original paper of this project), he should employ the CMDeo algorithm (instead of CMRules) available in our modification of the SPMF library "https://github.com/segura2010/SPMF". Instructions for this task have been gathered in the linked github project.


 


#	STEP - 5 	POST-PREPROCESSING

                                     
These lines are instructions for running the pmml.py python script, which transforms the (.txt) output rules file into a (.pmml) format annotated with gene and probe names. Note: the process should be repeated individually in each of the generated output rule files.


Required Input files:

	- The correspondence file that matches gene names with number ids. It was generated during the step 1 and it is available at "./GeneSeqRules/Output_files/Correspondence_file/Correspondence_gene_names"

	- The output rules files. They are available as: "./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/Output_rules_raw_group_1_measures.txt" for group 1 and "./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/Output_rules_raw_group_2_measures.txt" for group 2.


Instructions:

	(1) Open an Ubuntu-console or a Windows-Command Prompt.

	(2) Type the next lines: (for that purpose, you will need the latest python version installed)

	python ./GeneSeqRules/required_codes/pmml.py ./GeneSeqRules/Output_files/Correspondence_file/Correspondence_gene_names ./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/Output_rules_raw_group_1_measures.txt > ./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/pmml/Output_rules_raw_group_1_measures.pmml  

	for group 1 AND

	python ./GeneSeqRules/required_codes/pmml.py ./GeneSeqRules/Output_files/Correspondence_file/Correspondence_gene_names ./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/Output_rules_raw_group_2_measures.txt > ./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/pmml/Output_rules_raw_group_2_measures.pmml  

	for group 2



 

#	STEP - 6 	BIOLOGICAL VALIDATION


These lines are instructions for running the main.py python script, which annotates output rules with Biological Quality Measures. Note: the process should be repeated individually in each of the generated output file rules.

Required Input files:

	- The previously generated pmml output rule files. They are available as "./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/pmml/Output_rules_raw_group_1_measures.pmml" for group 1 and "./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/pmml/Output_rules_raw_group_2_measures.pmml" for group 2

	- The TRRUST database v2. It is available at "./GeneSeqRules/Data/TRRUST_database/TRRUST_rawdata.human.csv")

	- A file with by-probe annotations for GO and KEGG terms. This file should be generated at demand (depending on the Affymetrix array platform employed) using the R script available in "./GeneSeqRules/required_codes/Part_6_generating_annotation_file_with_go_and_kegg_terms_by_probe.r". Once generated, the required by-probe annotations file will be available in the directory "./GeneSeqRules/Data/ae_annot/final/ae.annots.csv"


Instructions:


	(1) Open an Ubuntu-console or a Windows-Command Prompt. Change the current directory to the directory "./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/pmml"

	(2) Then, type the next lines: (for that purpose, you will need the latest python version installed)

	python ./GeneSeqRules/required_codes/main.py ./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/pmml/Output_rules_raw_group_1_measures.pmml ./GeneSeqRules/Data/ae_annot/final/ae.annots.csv ./GeneSeqRules/Data/TRRUST_database/TRRUST_rawdata.human.csv

	for group 1 AND

	python ./GeneSeqRules/required_codes/main.py ./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/pmml/Output_rules_raw_group_2_measures.pmml ./GeneSeqRules/Data/ae_annot/ae.annots.csv ./GeneSeqRules/Data/TRRUST_database/TRRUST_rawdata.human.csv

	for group 2


	(3) The script will generate three output files in the current directory:

		- A (.pmml) output rule file in which extracted rules are annotated with the Biological quality metrics in the same directory "./GeneSeqRules/Output_files/Output_rules/Output_rules_more_measures/pmml" . This file will be the input file for the next step; the graphic visualization of results.

		- A (.xml) where specific GO and KEGG matches by rule are shown.

		- A (.xml) where specific Transcription Factors (TF) matches by rule are showed.





#	STEP - 7 	GRAPHICAL REPRESENTATION OF RESULTS


These lines are instructions for running the "./GeneSeqRules/required_codes/Part_7_NetworkGraph_gene_rules_customisable_plots.r" R script, which performs hierarchical edge bundling visualization of extracted rules. Note: the process should be repeated individually in each of the generated output file rules.


Required Script:

	- "./GeneSeqRules/required_codes/Part_7_NetworkGraph_gene_rules_customisable_plots.r"


Instructions:


	(1) Open an Ubuntu-console or a Windows-Command Prompt.

	(2) Run R environment 3.6.0 or higher.

	(3) Install the required packages listed in the script.

	(4) Run the script.


Scripts modifications: The scripts are adapted to the current example problem. If you change the localization of input required files or their name, you will need to modify the code. The code is also adapted for changing the colours and leading variables of the plot. These lines are highlighted in the code with the comment "# Customisable".





