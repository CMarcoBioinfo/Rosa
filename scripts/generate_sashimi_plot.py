import os,sys
import argparse
import random
import pandas as pd

#Create directory
def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory,exist_ok=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-directory', dest='directory', help='Path to ouput directory')
    parser.add_argument('-ggsashimi', dest='ggsashimi', help='full path to executable ggsashimi')
    parser.add_argument('-bam', dest='bam', help='full path to interest bam')
    parser.add_argument('-event_file', dest='event_file', help='full path to event interest')
    parser.add_argument('-extend_bp', dest='extend_bp', help='number of bp add to extremity. Default = 50', type=int, default=50)
    parser.add_argument('-MinThresholdNbReads', dest='MinThresholdNbReads', help='Minimum reads in junctions for plot. Default, half of the reads of the sample of interest junction', type=int, default=-1)
    parser.add_argument('-list_bam', dest='list_bam', help='list of bam')
    parser.add_argument('-nb_samples', dest='nb_samples', help='number of random sample in plot. Default = 4', type=int, default=4)
    parser.add_argument('-gtf', dest='gtf', help='full path to gtf of reference', default=-1)
    parser.add_argument('-color', dest='color', help='palette of color for plot', default=-1)
    args = parser.parse_args()

    #import options
    directory = args.directory
    ggsashimi = args.ggsashimi
    bam = args.bam
    event_file = args.event_file
    extend_bp = args.extend_bp
    MinThresholdNbReads = args.MinThresholdNbReads
    list_bam = args.list_bam.split() if args.list_bam else []
    nb_samples = args.nb_samples
    gtf = args.gtf
    color = args.color
    python = sys.executable

    create_directory_if_not_exists(directory)

    extension = event_file.rsplit(".", 1)[1]

    if extension == "xlsx":
        df = pd.read_excel(event_file, engine='openpyxl')
    else :
        df = pd.read_csv(event_file,sep = "\t")


    if df.empty:
        print(f"Aucuns évènements. Le fichier {event_file} est vide.")
        return

    extension_bam = ".sorted.bam"
    path_bam = bam.rsplit("/",1)[0] + "/"

    filtered_list = list()
    name_bam = ""
    for i in list_bam:
        string = i + extension_bam
        if bam.endswith(string):
            name_bam = i
        else:
            filtered_list.append(i)

    print("######################")
    print(f'#Generate sashimi plot for {name_bam}...')
    print("######################")
    
    for index, row in df.iterrows():
        chr = row["chr"]
        start = int(row["start"]) - extend_bp
        end = int(row["end"]) + extend_bp
        pdf = directory + "/" + row["Conca"] + ".pdf"

        if MinThresholdNbReads == -1:
            MinThresholdNbReads = int(row[[name_bam]].iloc[0] / 2)
        random_bam = random.sample(filtered_list, nb_samples)
        random_bam.insert(0, name_bam)

        bams = ""
        for j in random_bam :
            if j == name_bam :
                bams += j + "\t" + path_bam + j + extension_bam + "\tInterest" + "\n"
            else : 
                bams += j + "\t" + path_bam + j + extension_bam + "\tRandom" + "\n"
        bams = bams.rstrip("\n")

        cmd = f""" bash -c '{python} {ggsashimi} -b <(echo -e "{bams}") -c {chr}:{start}-{end} -M {MinThresholdNbReads} --fix-y-scale --height 2.75 --ann-height 2.75 --alpha 0.5 -o {pdf} """
        
        if gtf != -1:
            cmd += f"""-g {gtf} """ 

        if color != -1:
            cmd += f"""-P {color} -C 3 """
        
        cmd += "'"
        MinThresholdNbReads = args.MinThresholdNbReads
        os.system(cmd)

if __name__ == "__main__":
    main()


