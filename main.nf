#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { publish } from './modules/utils'
include { split } from './modules/utils'
include { concatenate } from './modules/utils'
include { translate2col } from './modules/utils'

include { dl_human } from './modules/uniprot_id_dict'
include { uniprotac_to_ensp_dict } from './modules/uniprot_id_dict'
include { IDa2uniprot2IDb } from './modules/uniprot_id_dict'

include { dl_string_aliases } from './modules/string_id_dict'

include { psp_kinase_substrates } from './modules/psp'
include { kinase_substrates_histograms } from './modules/psp'
include { kin_phos_clusters } from './modules/psp'
include { translate_kin_phos_clusters } from './modules/psp'
include { hs_phosphoproteome } from './modules/psp'

include { parse_hernandez2017_dataset } from './modules/hernandez2017'
include { translate_hernandez2017_metadata } from './modules/hernandez2017'
include { hernandez2017_rnk_translate1col_w_dict } from './modules/hernandez2017'
include { hernandez2017_make_seqrnk } from './modules/hernandez2017'
include { hernandez2017_make_uniprot_seqrnk } from './modules/hernandez2017'
include { hernandez2017_make_rnk } from './modules/hernandez2017'
include { benchmark_phosx_hernandez2017 } from './modules/hernandez2017'
include { benchmark_phosx_per_kinase_hernandez2017 } from './modules/hernandez2017'
include { link_kinase_annotations_to_metadata } from './modules/hernandez2017'
include { visualize_kinase_annotations } from './modules/hernandez2017'

include { parse_cptac_dataset } from './modules/cptac'
include { split_cptac_samples } from './modules/cptac'
include { make_cptac_rnk } from './modules/cptac'
include { make_cptac_seqrnk } from './modules/cptac'
include { make_cptac_uniprot_rnk } from './modules/cptac'
include { make_cptac_uniprot_seqrnk } from './modules/cptac'
include { translate_cptac_metadata } from './modules/cptac'
include { benchmark_phosx_cptac } from './modules/cptac'
include { benchmark_phosx_per_kinase_cptac } from './modules/cptac'

include { parse_kinomics_tremetinib_dataset } from './modules/kinomics'
include { parse_kinomics_tremetinib_metadata } from './modules/kinomics'
include { translate_kinomics_tremetinib_metadata } from './modules/kinomics'
include { parse_kinomics_vemurafenib_dataset } from './modules/kinomics'
include { parse_kinomics_vemurafenib_metadata } from './modules/kinomics'
include { translate_kinomics_vemurafenib_metadata } from './modules/kinomics'
include { join_kinomics_metadata } from './modules/kinomics'
include { benchmark_phosx_kinomics } from './modules/kinomics'
include { benchmark_phosx_per_kinase_kinomics } from './modules/kinomics'

include { run_phosx } from './modules/phosx'
include { run_phosx_serthr_only } from './modules/phosx'
include { run_phosx_nouae } from './modules/phosx'
include { run_phosx_serthr_only_nouae } from './modules/phosx'

include { get_kinex_scoring_matrix } from './modules/kinex'
include { run_kinex } from './modules/kinex'

include { make_kstar_networks } from './modules/kstar'
include { run_kstar } from './modules/kstar'
include { run_kstar_tolerant } from './modules/kstar'
include { translate_kstar_output } from './modules/kstar'

include { run_gsea } from './modules/gsea'

include { ptmsea_hernandez2017 } from './modules/ptmsea'
include { ptmsea_cptac } from './modules/ptmsea'
include { translate_ptmsea_output } from './modules/ptmsea'

include { zscore_hernandez2017 } from './modules/zscore'
include { zscore_cptac } from './modules/zscore'
include { translate_zscore_output } from './modules/zscore'


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


workflow GET_STRING_ID_DICT {

    main:
        string_id_dict = dl_string_aliases()

    emit:
        string_id_dict

}


/*workflow PSP {

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

}*/


workflow PSP {

    take:
        psp_kinase_substrates_dataset_gz
        psp_phosphosites_dataset_gz
        string_id_dict

    main:
        kinase_substrates = psp_kinase_substrates( psp_kinase_substrates_dataset_gz )

        // make histograms of kinase and phosphosites frequencies in kinase-substrate pairs
        kinase_substrates_histograms( kinase_substrates )

        kin_sub_clusters_untranslated = kin_phos_clusters( kinase_substrates )
        kin_sub_clusters = translate_kin_phos_clusters( kin_sub_clusters_untranslated,
                                                        string_id_dict )
        human_phosphosites = hs_phosphoproteome( psp_phosphosites_dataset_gz )

    emit:
        kin_sub_clusters
        human_phosphosites

}


workflow KINEX_SCORING_MATRIX {

    take:
        dict

    main:
        scoring_matrix = get_kinex_scoring_matrix( dict )

    emit:
        scoring_matrix

}


workflow HERNANDEZ2017_DATASET {

    take:
        uniprotac2ENSP_dict
        geneSynonym2geneName
        gene_synonym_2_gene_name_dict

    main:
        dataset = parse_hernandez2017_dataset()
        untranslated_rnk = dataset.rnk
                .flatMap()
                .map{file -> tuple( file.baseName, file )}
        //metadata = translate2col(dataset.metadata, geneSynonym2geneName )
        metadata = translate_hernandez2017_metadata( dataset.metadata, geneSynonym2geneName )
        metadata_annotated = link_kinase_annotations_to_metadata( metadata,
                                                                  gene_synonym_2_gene_name_dict,
                                                                  kinase_metadata_h5 )
        visualize_metadata = visualize_kinase_annotations( metadata_annotated )
        uniprot_rnk = hernandez2017_rnk_translate1col_w_dict( untranslated_rnk, uniprotac2ENSP_dict )
        seqrnk = hernandez2017_make_seqrnk( uniprot_rnk )
        uniprot_seqrnk = hernandez2017_make_uniprot_seqrnk( uniprot_rnk )
        rnk = hernandez2017_make_rnk( uniprot_rnk )
        publish( Channel.of('datasets/hernandez2017/metadata_translated.tsv').combine(metadata) )

    emit:
        metadata
        uniprot_rnk
        seqrnk
        uniprot_seqrnk
        rnk

}


workflow CPTAC_DATASET {

    take:
        uniprotac2ENSP_dict
        geneSynonym2geneName

    main:
        dataset = parse_cptac_dataset()
        data = split_cptac_samples( dataset.data )
                    .flatMap()
                    .map{file -> tuple( file.baseName, file )}
        rnk = make_cptac_rnk( data )
        uniprot_rnk = make_cptac_uniprot_rnk( rnk, uniprotac2ENSP_dict )
        seqrnk = make_cptac_seqrnk( data )
        uniprot_seqrnk = make_cptac_uniprot_seqrnk( uniprot_rnk )

        //metadata = translate2col(dataset.metadata, geneSynonym2geneName )
        metadata = translate_cptac_metadata( dataset.metadata, geneSynonym2geneName )
        publish( Channel.of('datasets/cptac/metadata_translated.tsv').combine(metadata) )

    emit:
        metadata
        uniprot_rnk
        seqrnk
        uniprot_seqrnk
        rnk

}


workflow KINOMICS_TREMETINIB_DATASET {

    take:
        geneSynonym2geneName

    main:
        dataset = parse_kinomics_tremetinib_dataset()
        seqrnk = dataset.seqrnk.flatMap()
                .map{file -> tuple( file.baseName, file )}
        uniprot_seqrnk = dataset.uniprot_seqrnk.flatMap()
                .map{file -> tuple( file.baseName, file )}
        rnk = dataset.rnk.flatMap()
                .map{file -> tuple( file.baseName, file )}
        metadata_untr = parse_kinomics_tremetinib_metadata()
        metadata = translate_kinomics_tremetinib_metadata( metadata_untr, geneSynonym2geneName )

    emit:
        seqrnk
        uniprot_seqrnk
        rnk
        metadata

}


workflow KINOMICS_VEMURAFENIB_DATASET {

    take:
        geneSynonym2geneName

    main:
        dataset = parse_kinomics_vemurafenib_dataset()
        seqrnk = dataset.seqrnk.flatMap()
                .map{file -> tuple( file.baseName, file )}
        uniprot_seqrnk = dataset.uniprot_seqrnk.flatMap()
                .map{file -> tuple( file.baseName, file )}
        rnk = dataset.rnk.flatMap()
                .map{file -> tuple( file.baseName, file )}
        metadata_untr = parse_kinomics_vemurafenib_metadata()
        metadata = translate_kinomics_vemurafenib_metadata( metadata_untr, geneSynonym2geneName )

    emit:
        seqrnk
        uniprot_seqrnk
        rnk
        metadata

}


workflow KINOMICS_DATASET {

    take:
        tremetinib_seqrnk
        tremetinib_uniprot_seqrnk
        tremetinib_rnk
        tremetinib_metadata
        vemurafenib_seqrnk
        vemurafenib_uniprot_seqrnk
        vemurafenib_rnk
        vemurafenib_metadata

    main:
        seqrnk = tremetinib_seqrnk.concat( vemurafenib_seqrnk )
        uniprot_seqrnk = tremetinib_uniprot_seqrnk.concat( vemurafenib_uniprot_seqrnk )
        rnk = tremetinib_rnk.concat( vemurafenib_rnk )
        metadata = join_kinomics_metadata( tremetinib_metadata, vemurafenib_metadata)

    emit:
        seqrnk
        uniprot_seqrnk
        rnk
        metadata

}


workflow PHOSX_HERNANDEZ2017 {

    take:
        seqrnk
        dict

    main:
        phosx_output = run_phosx( seqrnk, dict )
                            .tsv
                            .collect()

    emit:
        phosx_output

}


workflow PHOSXNOUAE_HERNANDEZ2017 {

    take:
        seqrnk
        dict

    main:
        phosx_output = run_phosx_nouae( seqrnk, dict )
                            .tsv
                            .collect()

    emit:
        phosx_output

}


workflow PHOSX_CPTAC {

    take:
        seqrnk
        dict

    main:
        phosx_output = run_phosx( seqrnk, dict )
                            .tsv
                            .collect()

    emit:
        phosx_output

}


workflow PHOSXNOUAE_CPTAC {

    take:
        seqrnk
        dict

    main:
        phosx_output = run_phosx_nouae( seqrnk, dict )
                            .tsv
                            .collect()

    emit:
        phosx_output

}


workflow PHOSX_KINOMICS {

    take:
        seqrnk
        dict

    main:
        phosx_output = run_phosx_serthr_only( seqrnk, dict )
                            .tsv
                            .collect()

    emit:
        phosx_output

}

workflow PHOSXNOUAE_KINOMICS {

    take:
        seqrnk
        dict

    main:
        phosx_output = run_phosx_serthr_only_nouae( seqrnk, dict )
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


workflow KINEX_CPTAC {

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


workflow KINEX_KINOMICS {

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


workflow GSEA_CPTAC {

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


workflow GSEA_KINOMICS {

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


workflow KSTAR_HERNANDEZ2017{

    take:
        string_id_dict
        uniprot_rnk

    main:
        kstar_networks = make_kstar_networks()
        kstar_output_untr = run_kstar( uniprot_rnk,
                                       kstar_networks.st_net,
                                       kstar_networks.y_net )
        kstar_output = translate_kstar_output( kstar_output_untr,
                                               string_id_dict )
                            .collect()

    emit:
        kstar_output

}


workflow KSTAR_CPTAC{

    take:
        string_id_dict
        uniprot_rnk

    main:
        kstar_networks = make_kstar_networks()
        kstar_output_untr = run_kstar( uniprot_rnk,
                                       kstar_networks.st_net,
                                       kstar_networks.y_net )
        kstar_output = translate_kstar_output( kstar_output_untr,
                                               string_id_dict )
                            .collect()

    emit:
        kstar_output

}


workflow KSTAR_KINOMICS{

    take:
        string_id_dict
        uniprot_rnk

    main:
        kstar_networks = make_kstar_networks()
        kstar_output_untr = run_kstar_tolerant( uniprot_rnk,
                                                kstar_networks.st_net,
                                                kstar_networks.y_net )
        kstar_output = translate_kstar_output( kstar_output_untr,
                                               string_id_dict )
                            .collect()

    emit:
        kstar_output

}


workflow PTMSEA_HERNANDEZ2017{

    take:
        string_id_dict

    main:
        ptmsea_output_untr = ptmsea_hernandez2017()
                                .flatMap()
                                .map{file -> tuple( file.baseName, file )}
        ptmsea_output = translate_ptmsea_output( ptmsea_output_untr,
                                                 string_id_dict )
                            .collect()

    emit:
        ptmsea_output

}


workflow PTMSEA_CPTAC{

    take:
        string_id_dict

    main:
        ptmsea_output_untr = ptmsea_cptac()
                                .flatMap()
                                .map{file -> tuple( file.baseName, file )}
        ptmsea_output = translate_ptmsea_output( ptmsea_output_untr,
                                                 string_id_dict )
                            .collect()

    emit:
        ptmsea_output

}


workflow ZSCORE_HERNANDEZ2017{

    take:
        string_id_dict

    main:
        zscore_output_untr = zscore_hernandez2017()
                                .flatMap()
                                .map{file -> tuple( file.baseName, file )}
        zscore_output = translate_zscore_output( zscore_output_untr,
                                                 string_id_dict )
                            .collect()

    emit:
        zscore_output

}


workflow ZSCORE_CPTAC{

    take:
        string_id_dict

    main:
        zscore_output_untr = zscore_cptac()
                                .flatMap()
                                .map{file -> tuple( file.baseName, file )}
        zscore_output = translate_zscore_output( zscore_output_untr,
                                                 string_id_dict )
                            .collect()

    emit:
        zscore_output

}


workflow BENCHMARK_PHOSX_HERNANDEZ2017 {

    take:
        phosx_output
        gsea_output
        kinex_output
        kstar_output
        ptmsea_output
        zscore_output
        phosxnouae_output
        metadata

    main:
        benchmark_phosx_hernandez2017(phosx_output,
                                      gsea_output,
                                      kinex_output,
                                      kstar_output,
                                      ptmsea_output,
                                      zscore_output,
                                      phosxnouae_output,
                                      metadata )
        benchmark_phosx_per_kinase_hernandez2017(phosx_output,
                                                 gsea_output,
                                                 kinex_output,
                                                 kstar_output,
                                                 ptmsea_output,
                                                 zscore_output,
                                                 phosxnouae_output,
                                                 metadata )

}


workflow BENCHMARK_PHOSX_CPTAC {

    take:
        phosx_output
        gsea_output
        kinex_output
        kstar_output
        ptmsea_output
        zscore_output
        phosxnouae_output
        metadata

    main:
        benchmark_phosx_cptac(phosx_output,
                              gsea_output,
                              kinex_output,
                              kstar_output,
                              ptmsea_output,
                              zscore_output,
                              phosxnouae_output,
                              metadata )
        benchmark_phosx_per_kinase_cptac(phosx_output,
                                         gsea_output,
                                         kinex_output,
                                         kstar_output,
                                         ptmsea_output,
                                         zscore_output,
                                         phosxnouae_output,
                                         metadata )

}


workflow BENCHMARK_PHOSX_KINOMICS {

    take:
        phosx_output
        gsea_output
        kinex_output
        kstar_output
        phosxnouae_output
        metadata

    main:
        benchmark_phosx_kinomics(phosx_output,
                                 gsea_output,
                                 kinex_output,
                                 kstar_output,
                                 phosxnouae_output,
                                 metadata )

        /*benchmark_phosx_per_kinase_kinomics(phosx_output,
                                            gsea_output,
                                            kinex_output,
                                            kstar_output,
                                            phosxnouae_output,
                                            metadata )*/

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


    // get alias table from STRING
    string_id_dict = GET_STRING_ID_DICT()


    // get human k-p associations from psp
    psp_data = PSP( psp_kinase_substrate_gz,
                    psp_phosphosites_dataset_gz,
                    string_id_dict )
    psp_kin_sub_clusters = psp_data.kin_sub_clusters
    //psp_human_phosphosites = psp_data.human_phosphosites


    // process Kinex scoring matrix
    scoring_matrix = KINEX_SCORING_MATRIX( string_id_dict )


    // get hernandez2017 dataset
    hernandez2017 = HERNANDEZ2017_DATASET( gene_id_dict.uniprotac2ENSP_dict,
                                           string_id_dict,
                                           gene_synonym_2_gene_name_dict )


    // get cptac dataset
    cptac = CPTAC_DATASET( gene_id_dict.uniprotac2ENSP_dict,
                           string_id_dict )


    // get kinomics dataset
    kinomics_tremetinib = KINOMICS_TREMETINIB_DATASET( string_id_dict )
    kinomics_vemurafenib = KINOMICS_VEMURAFENIB_DATASET( string_id_dict )
    kinomics = KINOMICS_DATASET( kinomics_tremetinib.seqrnk,
                                 kinomics_tremetinib.uniprot_seqrnk,
                                 kinomics_tremetinib.rnk,
                                 kinomics_tremetinib.metadata,
                                 kinomics_vemurafenib.seqrnk,
                                 kinomics_vemurafenib.uniprot_seqrnk,
                                 kinomics_vemurafenib.rnk,
                                 kinomics_vemurafenib.metadata )

    // run methods on hernandez2017
    phosx_hernandez2017 = PHOSX_HERNANDEZ2017( hernandez2017.seqrnk,
                                               string_id_dict )
    phosxnouae_hernandez2017 = PHOSXNOUAE_HERNANDEZ2017( hernandez2017.seqrnk,
                                                         string_id_dict )
    kinex_hernandez2017 = KINEX_HERNANDEZ2017( hernandez2017.seqrnk,
                                               string_id_dict,
                                               scoring_matrix )
    gsea_hernandez2017 = GSEA_HERNANDEZ2017( psp_kin_sub_clusters,
                                             hernandez2017.rnk )
    kstar_hernandez2017 = KSTAR_HERNANDEZ2017( string_id_dict,
                                               hernandez2017.uniprot_seqrnk )
    ptmsea_hernandez2017 = PTMSEA_HERNANDEZ2017( string_id_dict )
    zscore_hernandez2017 = ZSCORE_HERNANDEZ2017( string_id_dict )


    // run methods on cptac
    phosx_cptac = PHOSX_CPTAC( cptac.seqrnk,
                               string_id_dict )
    phosxnouae_cptac = PHOSXNOUAE_CPTAC( cptac.seqrnk,
                                         string_id_dict )
    kinex_cptac = KINEX_CPTAC( cptac.seqrnk,
                               string_id_dict,
                               scoring_matrix )
    gsea_cptac = GSEA_CPTAC( psp_kin_sub_clusters,
                             cptac.uniprot_rnk )
    kstar_cptac = KSTAR_CPTAC( string_id_dict,
                               cptac.uniprot_seqrnk )
    ptmsea_cptac = PTMSEA_CPTAC( string_id_dict )
    zscore_cptac = ZSCORE_CPTAC( string_id_dict )


    // run methods on Kinomics
    phosx_kinomics = PHOSX_KINOMICS( kinomics.seqrnk,
                                     string_id_dict )
    phosxnouae_kinomics = PHOSXNOUAE_KINOMICS( kinomics.seqrnk,
                                               string_id_dict )
    kinex_kinomics = KINEX_KINOMICS( kinomics.seqrnk,
                                     string_id_dict,
                                     scoring_matrix )
    gsea_kinomics = GSEA_KINOMICS( psp_kin_sub_clusters,
                                   kinomics.rnk )
    kstar_kinomics = KSTAR_KINOMICS( string_id_dict,
                                     kinomics.uniprot_seqrnk )
    

    // performance comparison on Hernandez2017
    BENCHMARK_PHOSX_HERNANDEZ2017( phosx_hernandez2017,
                                    gsea_hernandez2017,
                                    kinex_hernandez2017,
                                    kstar_hernandez2017,
                                    ptmsea_hernandez2017,
                                    zscore_hernandez2017,
                                    phosxnouae_hernandez2017,
                                    hernandez2017.metadata )


    // performance comparison on cptac
    BENCHMARK_PHOSX_CPTAC( phosx_cptac,
                            gsea_cptac,
                            kinex_cptac,
                            kstar_cptac,
                            ptmsea_cptac,
                            zscore_cptac,
                            phosxnouae_cptac,
                            cptac.metadata )


    // performance comparison on kinomics
    /*BENCHMARK_PHOSX_KINOMICS( phosx_kinomics,
                              gsea_kinomics,
                              kinex_kinomics,
                              kstar_kinomics,
                              phosxnouae_kinomics,
                              kinomics.metadata )*/

    // publish nextflow config
    PUBLISH_CONFIG()

}
