#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { publish } from './modules/utils'

include { dl_human } from './modules/uniprot_id_dict'
include { uniprotac_to_ensp_dict } from './modules/uniprot_id_dict'
include { IDa2uniprot2IDb } from './modules/uniprot_id_dict'

include { get_ser_thr_kinome_2023_suppl_table_2 } from './modules/pssm'
include { get_ser_thr_kinase_pssm } from './modules/pssm'

include { psp_kinase_substrates } from './modules/psp'
include { kin_phos_clusters } from './modules/psp'
include { translate_kin_phos_clusters } from './modules/psp'
include { hs_phosphoproteome } from './modules/psp'

include { parse_hernandez2017_dataset } from './modules/hernandez2017'
include { hernandez2017_rnk_translate1col_w_dict } from './modules/hernandez2017'
include { hernandez2017_make_seqrnk } from './modules/hernandez2017'
include { hernandez2017_make_rnk } from './modules/hernandez2017'
include { benchmark_phosx_hernandez2017 } from './modules/hernandez2017'

include { pssm_background_scores } from './modules/phossea'

include { run_phosx } from './modules/phosx'

include { dl_kinex_scoring_matrix } from './modules/kinex'
include { run_kinex } from './modules/kinex'

include { run_gsea } from './modules/gsea'

// ===== //

workflow GET_GENE_ID_DICT {

    main:
        uniprot_id_dict = dl_human()
        uniprotac2ENSP_dict = uniprotac_to_ensp_dict( uniprot_id_dict )

    emit:
        uniprotac2ENSP_dict
        uniprot_id_dict

}

workflow Gene_Synonym__2__Gene_Name {

    take:
        uniprot_id_dict

    main:
        dict = IDa2uniprot2IDb( uniprot_id_dict,
                                'Gene_Synonym',
                                'Gene_Name' )

    emit:
        dict

}

workflow SER_THR_KINASES_PSSM {

    take:
        gene_synonym_2_gene_name_dict

    main:
        kinome_2023_suppl_table_2 = get_ser_thr_kinome_2023_suppl_table_2()
        pssm_out = get_ser_thr_kinase_pssm( kinome_2023_suppl_table_2, gene_synonym_2_gene_name_dict )
        pssm = pssm_out.pssm_dict_h5

    emit:
        pssm

}

workflow PSP {

    take:
        psp_kinase_substrates_dataset_gz
        psp_phosphosites_dataset_gz
        dict

    main:
        kinase_substrates = psp_kinase_substrates( psp_kinase_substrates_dataset_gz )
        kin_sub_clusters_untranslated = kin_phos_clusters( kinase_substrates )
        kin_sub_clusters = translate_kin_phos_clusters( kin_sub_clusters_untranslated, dict )
        human_phosphosites = hs_phosphoproteome( psp_phosphosites_dataset_gz )

    emit:
        kinase_substrates
        kin_sub_clusters
        human_phosphosites

}

workflow PSSM_BACKGROUND_SCORES {

    take:
        human_phosphosites_tsv
        pssm_dict_h5

    main:
        pssm_bg_scores = pssm_background_scores( human_phosphosites_tsv,
                                                 pssm_dict_h5 )
        tsvgz = pssm_bg_scores.tsvgz
        h5 = pssm_bg_scores.h5

    emit:
        tsvgz
        h5

}

workflow KINEX_SCORING_MATRIX {

    take:
        dict

    main:
        scoring_matrix = dl_kinex_scoring_matrix(dict)

    emit:
        scoring_matrix

}

workflow HERNANDEZ2017_DATASET {

    take:
        uniprotac2ENSP_dict

    main:
        dataset = parse_hernandez2017_dataset()
        untranslated_rnk = dataset.rnk
                .flatMap()
                .map{file -> tuple( file.baseName, file )}
        metadata = dataset.metadata
        uniprot_rnk = hernandez2017_rnk_translate1col_w_dict( untranslated_rnk, uniprotac2ENSP_dict )
        seqrnk = hernandez2017_make_seqrnk( uniprot_rnk )
        rnk = hernandez2017_make_rnk( uniprot_rnk )

    emit:
        metadata
        uniprot_rnk
        seqrnk
        rnk

}

workflow PHOSX_HERNANDEZ2017 {

    take:
        seqrnk
        pssm_dict_h5
        pssm_bg_scores

    main:
        phosx_output = run_phosx( seqrnk, pssm_dict_h5, pssm_bg_scores )
                            .tsv
                            .collect()

    emit:
        phosx_output

}

workflow KINEX_HERNANDEZ2017 {

    take:
        seqrnk
        dict
        scoring_matrix
    
    main:
        kinex_output = run_kinex( seqrnk, dict, scoring_matrix )
                            .tsv
                            .collect()

    emit:
        kinex_output

}

workflow GSEA_HERNANDEZ2017 {

    take:
        psp_kin_sub_clusters
        rnk

    main:
        gsea_output = run_gsea( psp_kin_sub_clusters, rnk )
                        .csv
                        .collect()

    emit:
        gsea_output

}

workflow BENCHMARK_PHOSX_HERNANDEZ2017 {

    take:
        phosx_output
        gsea_output
        kinex_output
        metadata

    main:
        benchmark_phosx_hernandez2017(phosx_output,
                                      gsea_output,
                                      kinex_output,
                                      metadata )

}

workflow PUBLISH_CONFIG {

    main:
        config_ch = Channel.fromPath("${projectDir}/nextflow.config")
        val_ch = Channel.of('workflow/nextflow.config')
        
        publish( val_ch.combine( config_ch ) )

}

workflow {

    // get UniProt-based ID mapping dictionary: UniProtAc --> ENSP
    gene_id_dict = GET_GENE_ID_DICT()

    // get UniProt-based ID mapping dictionary: Gene Synonym --> Gene Name
    gene_synonym_2_gene_name_dict = Gene_Synonym__2__Gene_Name( gene_id_dict.uniprot_id_dict )

    // get pssms
    ser_thr_kinases_pssm_dict_h5 = SER_THR_KINASES_PSSM( gene_synonym_2_gene_name_dict ).pssm

    // get human k-p associations from psp
    psp_data = PSP( psp_kinase_substrate_gz,
                    psp_phosphosites_dataset_gz,
                    gene_synonym_2_gene_name_dict )
    psp_kin_sub_clusters = psp_data.kin_sub_clusters
    psp_human_phosphosites = psp_data.human_phosphosites

    // compute a priori pssm score distributions
    pssm_bg_scores = PSSM_BACKGROUND_SCORES( psp_human_phosphosites,
                                             ser_thr_kinases_pssm_dict_h5 )

    scoring_matrix = KINEX_SCORING_MATRIX( gene_synonym_2_gene_name_dict )
    
    // get hernandez2017 dataset
    hernandez2017 = HERNANDEZ2017_DATASET( gene_id_dict.uniprotac2ENSP_dict )
    
    // run methods on hernandez2017
    phosx_hernandez2017 = PHOSX_HERNANDEZ2017( hernandez2017.seqrnk,
                                               ser_thr_kinases_pssm_dict_h5,
                                               pssm_bg_scores.tsvgz )
    kinex_hernandez2017 = KINEX_HERNANDEZ2017( hernandez2017.seqrnk,
                                               gene_synonym_2_gene_name_dict,
                                               scoring_matrix )
    gsea_hernandez2017 = GSEA_HERNANDEZ2017( psp_kin_sub_clusters,
                                             hernandez2017.rnk )
    
    // performance comparison
    BENCHMARK_PHOSX_HERNANDEZ2017( phosx_hernandez2017,
                                   gsea_hernandez2017,
                                   kinex_hernandez2017,
                                   hernandez2017.metadata )


    PUBLISH_CONFIG()

}
