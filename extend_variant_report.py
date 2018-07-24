import requests
import json

#Variant class stores details from the input spreadsheet
class Variant:
    def __init__ (self, gene, cds_change, aa_change, preferred_transcript, hgvs_format):
        self.gene = gene
        self.cds_change = cds_change
        self.aa_change = aa_change
        self.preferred_transcript = preferred_transcript
        self.hgvs_format = hgvs_format

#removes formatting from protein_descriptions and transcript_descriptions
def clean_up (string):
    cut = string.split(':')
    half_clean = cut[1].replace("(", "")
    cleaned = half_clean.replace(")", "")
    return cleaned

def generate_output_file(infile):
    temp_string = infile.replace(".txt", "")
    return(temp_string+".MUTALYZER.txt")

#populates the preferred transcripts dictionary used to generate the query term for Mutalyzer
def read_preferred_transcripts(infile): 
    preferred_transcripts = {}
    with open(infile, 'r') as f:
        for line in f:
            fields = line.split("\t")
            preferred_transcripts[fields[0]] = fields[1]
    return (preferred_transcripts)


#creates a list of Variant objects that stores pertinent information from the input file
def create_variant_list(preferred_transcripts, infile):   
    #create a list to store all Variant objects
    variants = []
    
    #loop through the input file, instantiating and saving to the list Variant objects
    with open(infile, 'r') as f:
        for line in f:
            if (line.startswith("Chr")): #skip the header row
                continue
            fields = line.split("\t") #list of each column
            values = [] #temp array values stores the information as it is being read to later instantiate Variant objects
            gene = fields[2]
            values.append(gene) #gene 
            cds_change = fields[13]
            values.append(cds_change) #cds_change
            values.append(fields[15]) #aa_change
            try:
                preferred_transcript = preferred_transcripts[gene].rstrip() #rstrip removes trailing whitespace for proper formatting
                values.append(preferred_transcript)
            except KeyError:
                values.append("NOT LOCATED")
            hgvs_format = preferred_transcript + ":c." + cds_change
            values.append(hgvs_format)
            item = Variant(*values) #Instantiates Variant items
            variants.append(item) #Adds Tile tiem to the list "variants"              
    return variants

 #populates lists that represent each column (with header) containing each set of values to add
def create_lists(variants):
    url1 = "https://mutalyzer.nl/json/runMutalyzerLight?variant="
    url2 = "https://mutalyzer.nl/json/numberConversion?build=hg19;variant="
    
    #create lists containing all values for each column to add to the spreadsheet
    protein_descriptions = []
    transcript_descriptions = []
    version_numbers = []   
    chromosomal_locations = []
    
    #Set the first value of each column to the header value
    protein_descriptions.append("AA Change (Mutalyzer)")
    transcript_descriptions.append("CDS Change (Mutalyzer)")
    version_numbers.append("Ref Transcript.ncbi")
    chromosomal_locations.append("Chromosomal Location (Mutalyzer)")
    
    #populate lists containing all values in each column to add to the spreadsheet
    idx = 0
    for i in variants:
        #ping Mutalyzer's runMutalyzerLight's function
        ping1_url = url1 + variants[idx].hgvs_format
        r = requests.get(ping1_url, verify = False)
       
        #save into a dictionary called dict
        dict = json.loads(r.text)
    
        #parse individual elements
        if (dict['proteinDescriptions'] == []):
            protein_description = ""
        else:
            protein_description = clean_up(dict['proteinDescriptions'][0]) 
        
        if (dict['transcriptDescriptions'] == []):
            transcript_description = ""
        else:
            transcript_description = clean_up(dict['transcriptDescriptions'][0]) 
        
        if(dict['legend'] == []):
            version_number = ""
        else:
           version_number = dict['legend'][0]['id']
        
        #populate lists
        protein_descriptions.append(protein_description)
        transcript_descriptions.append(transcript_description)
        version_numbers.append(version_number)
        
        #ping Mutalyzer's numberConversion function
        ping2_url = url2 + version_number + ":" + transcript_description
        r = requests.get(ping2_url, verify = False)
        result = json.loads(r.text)
        
        #parse element to save g.
        if(result[0] == None):
           chromosomal_location = ""
        elif(result[0] == ""):
           chromosomal_location = ""
        else:
           split = result[0].split(":")
           chromosomal_location = split[1]
        
        #populate list
        chromosomal_locations.append(chromosomal_location)
        
        #increment counter
        idx += 1       
    
    #save a master_list that contains all 4 new column lists to pass into the print_to_file function    
    master_list = [version_numbers, protein_descriptions, transcript_descriptions, chromosomal_locations]
        
    return (master_list)


#read to a file with the same name as the previous file + .MUTALYZER.txt
def print_to_file(read_file, master_list):
    #call to external function to create output file
    output_file = generate_output_file(read_file)
    
    #loop through reading file to append 4 new columns to each line
    idx = 0
    with open (output_file, 'w') as wf:
        with open (read_file, 'r') as rf:
            for line in rf:
                    wf.write(line.rstrip() + '\t' + master_list[0][idx] + '\t' + master_list[1][idx] + '\t' + master_list[2][idx] + '\t' + master_list[3][idx] + '\n')
                    idx += 1       
        
#main routine needs file to read from and file from which to populate the preferred transcripts dictionary    
def main(preferred_transcripts_file, variants_file):
    preferred_transcripts = read_preferred_transcripts(preferred_transcripts_file)
    variants = create_variant_list(preferred_transcripts, variants_file)
    masterList = create_lists(variants)
    print_to_file(variants_file, masterList)
    
if __name__ == '__main__':
    main("N:\Clinlab\Private\Molecular Pathology\mpdocs\MDX Staff Files\Laila\heme_mutalyzer/Heme-STAMP_preferred_transcripts_180522.txt", "N:\Clinlab\Private\Molecular Pathology\mpdocs\MDX Staff Files\Laila\heme_mutalyzer/HD701_HEME0030.variant_report.txt")    



