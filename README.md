# Public Health Bioinformatics (PHB) - WAPHL
Bioinformatics workflows for characterization, epidemiology and sharing of pathogen genomes. This is the Washington State Department of Health Public Health Lab fork of the Public Health Bioinformatics repository created by Theiagen.

**More information about the steps undertaken in these workflows is available on the Wiki for this repository for all pipelines not containing WAPHL in the workflow filename or via the [Theiagen Public Resources Documentation](https://theiagen.notion.site/Theiagen-Public-Health-Resources-a4bd134b0c5c4fe39870e21029a30566) for all other pipelines.**

Support for running these workflows can be sought by raising a [GitHub issue](https://github.com/theiagen/public_health_bioinformatics/issues/new) or by contacting Theiagen at support@theiagen.com. For workflows ending with WAPHL please contact [Holly Halstead](mailto:holly.halstead@doh.wa.gov).

These workflows are written in [WDL](https://github.com/openwdl/wdl), a language for specifying data processing workflows with a human-readable and writeable syntax. They have been developed by [Theiagen Genomics](https://theiagen.com/) and WAPHL to primarily run on the [Terra.bio](https://terra.bio/) platform.

### Contributors & Influence
* Based on collaborative work with Andrew Lang, PhD & his [Genomic Analysis WDL workflows](https://github.com/AndrewLangvt/genomic_analyses)
* Workflows and task development influenced by The Broad's [Viral Pipes](https://github.com/broadinstitute/viral-pipelines)
* TheiaCoV workflows for viral genomic characterization influenced by UPHL's [Cecret](https://github.com/UPHL-BioNGS/Cecret) & StaPH-B's [Monroe](https://staph-b.github.io/staphb_toolkit/workflow_docs/monroe/)
* TheiaProk workflows for bacterial genomic characterization influenced by Robert Petit's [bactopia](https://github.com/bactopia/bactopia)
* The PHB workflow user community. To provide feedback to pipeline's not ending in WAPHL, please visit the original repository for [Public Health Bioinformatics](https://github.com/theiagen/public_health_bioinformatics).
* Tasks using in WAPHL pipelines are highly derived from tasks in the original [Public Health Bioinformatics](https://github.com/theiagen/public_health_bioinformatics) repository created by Theiagen


*“DISCLAIMER: This test was adopted, and its performance characteristics determined by the WA PHL. It has not been approved by the FDA. Results are for infection control and epidemiological purposes and should not be used as the sole means for clinical diagnosis or patient management.”*

### Citation

Please cite this paper if publishing work using any pipelines **not containing WAPHL in the workflow filename**:

Libuit, Kevin G., Emma L. Doughty, James R. Otieno, Frank Ambrosio, Curtis J. Kapsak, Emily A. Smith, Sage M. Wright, et al. 2023. “Accelerating Bioinformatics Implementation in Public Health.” Microbial Genomics 9 (7). https://doi.org/10.1099/mgen.0.001051.

Alternatively, please cite this paper if using the TheiaEuk workflow:

Ambrosio, Frank, Michelle Scribner, Sage Wright, James Otieno, Emma Doughty, Andrew Gorzalski, Danielle Siao, et al. 2023. “TheiaEuk: A Species-Agnostic Bioinformatics Workflow for Fungal Genomic Characterization.” Frontiers in Public Health 11. https://doi.org/10.3389/fpubh.2023.1198213.
