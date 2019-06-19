# Run Evoquer on a list of intervals.
#
# Description of inputs:
#
#  Required:
#    String gatk_docker                 - GATK Docker image in which to run
#
#    File ref_fasta_dict                - HG38 Reference sequence dictionary.
#
#    Array[String] intervals            - Variant Context File (VCF) containing the variants to annotate.
#    String output_file_base_name       - Base name of desired output file WITHOUT extension.
#
#  Optional:
#
#     Boolean compress (default: false) - if true, will compress resulting VCF output.
#
#     File gatk4_jar_override           -  Override Jar file containing GATK 4.  Use this when overriding the docker JAR or when using a backend without docker.
#     Int  mem_gb                       -  Amount of memory to give to the machine running each task in this workflow (in gb).
#     Int  preemptible_attempts         -  Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb                -  Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                          -  Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb            -  Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
workflow Evoquer {
    String gatk_docker

    File ref_dict

    Array[String] intervals
    String output_file_base_name

    Boolean compress = false

    File? gatk4_jar_override
    Int?  mem_gb
    Int?  preemptible_attempts
    Int?  disk_space_gb
    Int?  cpu
    Int?  boot_disk_size_gb

    call EvoquerTask {
        input:
            ref_dict              = ref_dict,
            intervals             = intervals,
            output_file_base_name = output_file_base_name,
            compress              = compress,

            gatk_docker           = gatk_docker,
            gatk_override         = gatk4_jar_override,
            mem                   = mem_gb,
            preemptible_attempts  = preemptible_attempts,
            disk_space_gb         = disk_space_gb,
            cpu                   = cpu,
            boot_disk_size_gb     = boot_disk_size_gb
    }

    output {
        File genotyped_variant_calls_file = EvoquerTask.evoked_variants
        File genotyped_variant_calls_file_idx = EvoquerTask.evoked_variants_index
    }
}

################################################################################

task EvoquerTask {

    # ------------------------------------------------
    # Input args:
    File ref_dict

    Array[String] intervals
    String output_file_base_name

    Boolean compress

    # Runtime Options:
    String gatk_docker
    File? gatk_override
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    # ------------------------------------------------
    # Misc settings:
    String dollar = "$"

    # ------------------------------------------------
    # Process input args:

    # Output Names:
    String output_file = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
    String output_file_index = output_file +  if compress then ".tbi" else ".idx"

    # Timing info:
    String timing_output_file = "Evoquer.timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 1024 * 3
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ${timing_output_file}

        # Run Evoquer:
        gatk --java-options "-Xmx${command_mem}m" \
            Evoquer \
                --intervals ${default="" sep=" --intervals " intervals} \
                --sequence-dictionary ${ref_dict} \
                -O ${output_file}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ${timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: select_first([cpu, 1])
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File evoked_variants = "${output_file}"
        File evoked_variants_index = "${output_file_index}"
    }
 }
