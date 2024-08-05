# 1997  featureCounts -p -T 2 -a ../../../../06-data/05-GTF/GRCh37.P13/GRCh37.p13.SLA.chr.gtf -o test.txt 202304-1409241167-SLA.sorted.bam
#Include script
from scripts import merge_counts as mc

rule featureCounts:
    input:
        bai = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai",
        annotation = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["DATA_INPUT"]["ANNOTATION"]

    output:
        summaryFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts.summary",
        countsFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts",
        metadata = temp(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/sample_names.txt")



    params:
        directory = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/",
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
        featureCounts = config["DEPENDANCES"]["FEATURECOUNTS"],
        reads = "{reads}"

    threads:
        config["PARAMS"]["FEATURECOUNTS"]["THREADS"]

    run:
        create_directory_if_not_exists(params["directory"])
        shell("{params.featureCounts} -p "
        "-T {threads} "
        "-a {input.annotation} "
        "-o {params.directory}{params.reads}.counts {params.bam}")
        lock_file(output["metadata"])
        add_file(params["reads"], output["countsFile"])
        unclock_file(output["metadata"])

    
rule finalize_count:
    input:
        countsFile = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts",reads=all_samples)
    
    output:
        file = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/.done.txt"

    shell:
        "touch {output.file}"


 

rule concatene_new_files:
    input:
        get_available_files()

    output:
        concatene = temp(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/.concatene.txt")
    
    run:
        if (os.path.exists(output["concatene"])):
            df_concatene = pd.read_csv(output["concatene"],sep="\t",index_col=0,comment="#")
        else:
            df_concatene = pd.DataFrame()

        for i, file in enumerate(input):
            if (os.path.exists(output["concatene"])) or i > 0:
                df = pd.read_csv(file, sep="\t", index_col=0, usecols=[0, 6], comment="#")
            else:
                df = pd.read_csv(file, sep="\t", index_col=0, comment="#")
            df_concatene = pd.concat([df_concatene,df], ignore_index=True)

        df_concatene.to_csv(output["concatene"], sep="\t", index=False)



rule monitor_and_merge:
    input:
        metadata = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/sample_names.txt"

    output:
        mergeFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/merge.counts"
    
    params:
        done = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/.done.txt",
        concatene = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/.concatene.txt"


    run:
        inotify = INotify()
        watch_flags = flags.MODIFY
        wd = inotify.add_watch(input["metadata"], watch_flags)
        delay = 5
        first_event_time = None
        
        while True:
            events = inotify.read(timeout = 60000)
            current_time = time.time()
            if events:
                for event in events:
                    if event.mask & flags.MODIFY:
                        if first_event_time is None:
                            first_event_time = current_time
                        print("Modification détectée,atttente pour regrouper les événements...")

            if first_event_time and current_time - first_event_time >= delay:
                print("Exécution de concatene_new_files après regroupement des événements...")
                os.system("snakemake --cores 1 concatene_new_files")
                first_event_time = None

            if os.path.exists(params["done"]):
                print("Tous les fichiers de comptage sont disponibles. Regroupement final...")
                os.system("snakemake --cores 1 concatene_new_files")
                shell("mv {params.concatene} {output.mergeFile} && rm {params.done}")
                break

