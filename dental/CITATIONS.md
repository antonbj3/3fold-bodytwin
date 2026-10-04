# Datasets and citations

This register covers external datasets used or named by the supplied source implementations and the published numerical tables included in the portable profiles. Anatomical datasets, patient-derived geometry, labels, records and model weights are **not included** in this repository. The default 23 profiles use synthetic inputs or the small published group summaries identified below; they do not download these datasets.

Upstream citation and licence metadata was checked on **4 October 2026**. Required Ditto citations below reproduce the page's APA citation field verbatim, including its spelling and bibliographic abbreviations. DOI links missing from that field are supplied separately from its BibTeX metadata or the primary publication. A citation is attribution, not permission to redistribute data. Apache-2.0 covers this repository's software; each external dataset retains its own terms.

Creative Commons obligations: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) requires attribution, a licence link and disclosure of changes. [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) additionally requires distributed adaptations to use the same or a compatible licence. [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) restricts use to noncommercial purposes and adds ShareAlike. [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/) restricts use to noncommercial purposes and does not permit sharing adapted material. Use the linked terms for the exact conditions and exceptions.

## Bits2Bites

**Repository use:** Registered upper/lower intraoral surfaces and occlusion labels for contact, crown and occlusion source pipelines (for example X2, X54 and GENCAD).

**Source:** [Dataset page](https://ditto.ing.unimore.it/bits2bites/).

**Licence and access:** CC BY-NC-SA 4.0, explicitly linked by the dataset page. Attribution, noncommercial use and ShareAlike for distributed adaptations. Account registration is required for download. No Bits2Bites data is bundled.

**Required citations, verbatim from the dataset page:**

> Borghi, L., Lumetti, L., Cremonini, F., Rizzo, F., Grana, C., Lombardo, L., & Bolelli, F. (2025). Bits2Bites: Intra-oral Scans Occlusal Classification. In Oral and Dental Image Analysis Workshop - MICCAI 2025.

**Published chapter DOI:** [10.1007/978-3-032-20711-1_5](https://doi.org/10.1007/978-3-032-20711-1_5). The later proceedings record does not replace the dataset page's requested 2025 citation.

## Bite2Text

**Repository use:** Paired intraoral surfaces and clinician report categories in the X7/X21 contact-report pipelines; the FALT_TANDLAST sources retain the absence of a matched Bits2Bites case linkage.

**Source:** [Dataset page](https://ditto.ing.unimore.it/bite2text/).

**Licence and access:** CC BY-NC-SA 4.0, explicitly linked by the current dataset page. This updates earlier source notes that lacked a verified licence; it does not retrospectively identify a particular local archive version. Account registration is required. Photographs, reports, scans and trained weights are not bundled. All three citations requested by the page are retained below.

**Required citations, verbatim from the dataset page:**

> Lumetti, L., Rizzo, F., Cremonini, F., Candeloro, E., Luca, L., Grana, C., & Bolelli, F. (2026). Do Multimodal LLMs Understand Intraoral Dental Data? Dataset, Platform, and Baselines. European Conference on Computer Vision (ECCV).

> Borghi, L., Lumetti, L., Cremonini, F., Rizzo, F., Grana, C., Lombardo, L., & Bolelli, F. (2025). Bits2Bites: Intra-oral Scans Occlusal Classification. In Oral and Dental Image Analysis Workshop - MICCAI 2025.

**Published chapter DOI:** [10.1007/978-3-032-20711-1_5](https://doi.org/10.1007/978-3-032-20711-1_5). The later proceedings record does not replace the dataset page's requested 2025 citation.

> Bolelli, F., Lumetti, L., Nistelrooij, N. V., Vinayahalingam, S., Bartolomeo, M. D., Marchesini, K., Pellacani, A., Candeloro, E., Rosati, G., Xi, T., Isensee, F., Kirchhoff, Y., Kraemer, L., Rokuss, M., Ulrich, C., Maier-Hein, K., Jiang, Y., Liu, Y., Wang, L., Wang, H., Chen, S., Cui, Z., Shi, P., Pan, Z., Liang, X., Ma, Q., Konukoglu, E., Wodzinski, M., Müller, H., Mai, H., Dang, X., Bhandary, S., Grosu, R., Bergé, S., Anesi, A., & Grana, C. (2026). Multi-structure segmentation in CBCT volumes: The ToothFairy2 challenge. Medical Image Analysis, 104095. https://doi.org/10.1016/j.media.2026.104095

**DOIs:** [`10.1016/j.media.2026.104095`](https://doi.org/10.1016/j.media.2026.104095).

The page's citation field omits its ECCV DOI; the [publisher chapter](https://link.springer.com/chapter/10.1007/978-3-032-37574-2_26) supplies **10.1007/978-3-032-37574-2_26**. The Bits2Bites chapter DOI is **10.1007/978-3-032-20711-1_5**. Its APA field spells one author “Luca, L.”; that spelling is preserved rather than silently corrected.

## ToothFairy (MICCAI 2023)

**Repository use:** Inferior alveolar canal segmentation and canal-distance/source-comparison pipelines, distinct from the larger ToothFairy2 label set.

**Source:** [Dataset page](https://ditto.ing.unimore.it/toothfairy/).

**Licence and access:** The [official challenge dataset page](https://toothfairy.grand-challenge.org/dataset/) states CC BY-SA; that page does not specify the version number. Preserve the actual archive terms rather than substituting an article licence. Account registration is required. The challenge test set is private and will not be released. No training or test scans are bundled. Publications using the data must explicitly reference the ToothFairy challenge and cite all three works below.

**Required citations, verbatim from the dataset page:**

> Bolelli, F., Lumetti, L., Vinayahalingam, S., Di Bartolomeo, M., Pellacani, A., Marchesini, K., ... & Grana, C. (2024). Segmenting the Inferior Alveolar Canal in CBCTs Volumes: the ToothFairy Challenge. IEEE Transactions on Medical Imaging.

> Lumetti, L., Pipoli, V., Bolelli, F., Ficarra, E., & Grana, C. (2024). Enhancing Patch-Based Learning for the Segmentation of the Mandibular Canal. IEEE Access, 1–12. https://doi.org/10.1109/ACCESS.2024.3408629

> Cipriano, M., Allegretti, S., Bolelli, F., Pollastri, F., & Grana, C. (2022). Improving Segmentation of the Inferior Alveolar Nerve through Deep Label Propagation. In IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR) (pp. 21105–21114). IEEE.

**DOIs:** [`10.1109/ACCESS.2024.3408629`](https://doi.org/10.1109/ACCESS.2024.3408629), [`10.1109/CVPR52688.2022.02046`](https://doi.org/10.1109/CVPR52688.2022.02046), [`10.1109/TMI.2024.3523096`](https://doi.org/10.1109/TMI.2024.3523096).

## ToothFairy2 (MICCAI 2024)

**Repository use:** CBCT labels for tooth, jaw, canal, implant and restoration geometry in source-only preparation, pulp-port and implant-revision operations (for example X95 and X98).

**Source:** [Dataset page](https://ditto.ing.unimore.it/toothfairy2/).

**Licence and access:** CC BY-SA 4.0 in the [official benchmark dataset.json](https://github.com/AImageLab-zip/ToothFairy2-Benchmark/blob/main/dataset.json); the [challenge dataset page](https://toothfairy2.grand-challenge.org/dataset/) confirms CC BY-SA. Account registration is required. The hidden challenge test set remains private. No ToothFairy2 scans or labels are bundled. Publications must reference the ToothFairy2 challenge and cite all three works below.

**Required citations, verbatim from the dataset page:**

> Bolelli, F., Lumetti, L., Nistelrooij, N. V., Vinayahalingam, S., Bartolomeo, M. D., Marchesini, K., Pellacani, A., Candeloro, E., Rosati, G., Xi, T., Isensee, F., Kirchhoff, Y., Kraemer, L., Rokuss, M., Ulrich, C., Maier-Hein, K., Jiang, Y., Liu, Y., Wang, L., Wang, H., Chen, S., Cui, Z., Shi, P., Pan, Z., Liang, X., Ma, Q., Konukoglu, E., Wodzinski, M., Müller, H., Mai, H., Dang, X., Bhandary, S., Grosu, R., Bergé, S., Anesi, A., & Grana, C. (2026). Multi-structure segmentation in CBCT volumes: The ToothFairy2 challenge. Medical Image Analysis, 104095. https://doi.org/10.1016/j.media.2026.104095

> Bolelli, F., Marchesini, K., van Nistelrooij, N., Lumetti, L., Pipoli, V., Ficarra, E., Vinayahalingam, S., & Grana, C. (2025). Segmenting Maxillofacial Structures in CBCT Volume. In IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR) (pp. 1–10). IEEE. https://dx.doi.org/10.1109/CVPR52734.2025.00494

> Lumetti, L., Pipoli, V., Bolelli, F., Ficarra, E., & Grana, C. (2024). Enhancing Patch-Based Learning for the Segmentation of the Mandibular Canal. IEEE Access, 1–12. https://doi.org/10.1109/ACCESS.2024.3408629

**DOIs:** [`10.1016/j.media.2026.104095`](https://doi.org/10.1016/j.media.2026.104095), [`10.1109/ACCESS.2024.3408629`](https://doi.org/10.1109/ACCESS.2024.3408629), [`10.1109/CVPR52734.2025.00494`](https://doi.org/10.1109/CVPR52734.2025.00494).

## Maxillo / IAN alveolar canal

**Repository use:** Earlier mandibular-canal benchmark/segmentation provenance and historical nerve-geometry input, named by the extracted source implementations.

**Source:** [Dataset page](https://ditto.ing.unimore.it/maxillo/).

**Licence and access:** UNKNOWN: the public dataset page requires citation and registration but gives no explicit data licence in the retrieved page. The code/article licence does not establish the data terms. Obtain the applicable download agreement before redistribution. Account registration is required. Maxillo is the older dataset superseded by ToothFairy. No CBCT volume, annotation or model checkpoint is bundled.

**Required citations, verbatim from the dataset page:**

> Cipriano, M., Allegretti, S., Bolelli, F., Di Bartolomeo, M., Pollastri, F., Pellacani, A., Minafra, P., Anesi, A., & Grana, C. (2022). Deep Segmentation of the Mandibular Canal: a New 3D Annotated Dataset of CBCT Volumes. IEEE Access, 10, 11500–11510.

> Cipriano, M., Allegretti, S., Bolelli, F., Pollastri, F., & Grana, C. (2022). Improving Segmentation of the Inferior Alveolar Nerve through Deep Label Propagation. In IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR) (pp. 21105–21114). IEEE.

**DOIs:** [`10.1109/ACCESS.2022.3144840`](https://doi.org/10.1109/ACCESS.2022.3144840), [`10.1109/CVPR52688.2022.02046`](https://doi.org/10.1109/CVPR52688.2022.02046).

## Pulpy3D

**Repository use:** Pulp-cavity and root-canal annotation inputs to preparation-shell, remaining-tissue and root-canal source pipelines; these are annotation-derived geometries rather than calibrated tissue measurements.

**Source:** [Dataset page](https://ditto.ing.unimore.it/pulpy3d/).

**Licence and access:** UNKNOWN: neither the retrieved dataset page nor the [authors' repository README](https://github.com/KnightsLab/Pulpy3D) provides an explicit data licence. A licence for the underlying ToothFairy images does not establish the terms for the additional pulp annotations. Obtain the applicable download agreement before redistribution. Ditto requires account registration. No Pulpy3D scans, pulp/canal labels or derived anatomical surfaces are bundled.

**Required citations, verbatim from the dataset page:**

> Gamal, M., Baraka, M., & Torki, M. (2024). Automatic Mandibular Semantic Segmentation of Teeth Pulp Cavity and Root Canals, and Inferior Alveolar Nerve on Pulpy3D Dataset. In Medical Image Computing and Computer Assisted Intervention -- MICCAI 2024 (pp. 14–23). Cham: Springer Nature Switzerland.

> Lumetti, L., Pipoli, V., Bolelli, F., Ficarra, E., & Grana, C. (2024). Enhancing Patch-Based Learning for the Segmentation of the Mandibular Canal. IEEE Access, 1–12. https://doi.org/10.1109/ACCESS.2024.3408629

**DOIs:** [`10.1007/978-3-031-72111-3_2`](https://doi.org/10.1007/978-3-031-72111-3_2), [`10.1109/ACCESS.2024.3408629`](https://doi.org/10.1109/ACCESS.2024.3408629).

## Teeth3DS / Teeth3DS+

**Repository use:** Intraoral tooth instance/FDI segmentation and the training provenance of tooth-labelled source pipelines (for example X11 and X21). No trained weights or dataset meshes are bundled.

**Source:** [Authors' dataset site](https://crns-smartvision.github.io/teeth3ds/) and [official challenge repository, License section](https://github.com/abenhamadou/3DTeethSeg_MICCAI_Challenges#license). The latter explicitly assigns **CC BY-NC-ND 4.0 to the data**: attribution, noncommercial use, and no sharing of adapted data. The dataset site's CC BY-SA 4.0 footer concerns the website and must not be substituted for the data licence. Keep the exact version's terms when acquiring the OSF release; no authentication requirement was established from the cited public README.

**Required citation block, verbatim from the authors' dataset site:**

```bibtex
@article{ben2022teeth3ds,
title={{Teeth3Ds+: An Extended Benchmark for Intra-oral 3D Scans Analysis}},
author={Ben-Hamadou, Achraf and Neifar, Nour and Rekik, Ahmed and Smaoui, Oussama and Bouzguenda, Firas and Pujades, Sergi and  Boyer, Edmond and Ladroit, Edouard},
journal={arXiv preprint arXiv:2210.06094},
year={2022}
}

@article{ben20233dteethseg,
title={3DTeethSeg'22: 3D Teeth Scan Segmentation and Labeling Challenge},
author={Achraf Ben-Hamadou and Oussama Smaoui and Ahmed Rekik and Sergi Pujades and Edmond Boyer and Hoyeon Lim and Minchang Kim and Minkyung Lee and Minyoung Chung and Yeong-Gil Shin and Mathieu Leclercq and Lucia Cevidanes and Juan Carlos Prieto and Shaojie Zhuang and Guangshun Wei and Zhiming Cui and Yuanfeng Zhou and Tudor Dascalu and Bulat Ibragimov and Tae-Hoon Yong and Hong-Gi Ahn and Wan Kim and Jae-Hwan Han and Byungsun Choi and Niels van Nistelrooij and Steven Kempers and Shankeeth Vinayahalingam and Julien Strippoli and Aurélien Thollot and Hugo Setbon and Cyril Trosset and Edouard Ladroit},
journal={arXiv preprint arXiv:2305.18277},
year={2023}
}
```

**Preprint DOIs:** [10.48550/arXiv.2210.06094](https://doi.org/10.48550/arXiv.2210.06094) and [10.48550/arXiv.2305.18277](https://doi.org/10.48550/arXiv.2305.18277). These do not change the data licence or establish the licence of a locally selected snapshot.

## STS-3D / STS-Tooth3D / STS-Tooth

**Repository use:** CBCT tooth-mask and imaging-resolution source experiments (for example X1B, X6 and X30); the dataset record covers both the STS-2D-Tooth and STS-3D-Tooth subsets.

**Source:** [Zenodo record](https://zenodo.org/records/10597292).

**Licence and access:** CC BY 4.0 in the record metadata: attribution, licence link and disclosure of changes. The files are marked open; no account registration requirement is stated. No image, mask, crown mesh or case record is bundled.

**Exact citation from the record's BibTeX export:**

```bibtex
@dataset{wang_2024_10597292,
  author       = {wang, yaqi},
  title        = {STS-Tooth: A multi-modal dental dataset for semi-
                   supervised deep learning image segmentation
                  },
  month        = sep,
  year         = 2024,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.10597292},
  url          = {https://doi.org/10.5281/zenodo.10597292},
}
```

**Dataset DOI:** [10.5281/zenodo.10597292](https://doi.org/10.5281/zenodo.10597292).

## AlignerMovement (Zenodo 11280343)

**Repository use:** X37/X48 planned-versus-observed registered crown-motion sequences; observations of two cases are not a population force calibration.

**Source:** [Zenodo record](https://zenodo.org/records/11280343).

**Licence and access:** CC BY 4.0 in the record metadata: attribution, licence link and disclosure of changes. The files are marked open; no account registration requirement is stated. No image, mask, crown mesh or case record is bundled.

**Exact citation from the record's BibTeX export:**

```bibtex
@dataset{tanner_2024_11280343,
  author       = {Tanner, Christine and
                  Filippon, Ignacio and
                  von Jackowski, Jeannette A. and
                  Schulz, Georg and
                  Toepper, Tino and
                  Müller, Bert},
  title        = {Determining aligner-induced tooth movements in
                   three dimensions using clinical data of two
                   patients: datasets
                  },
  month        = may,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v1.0},
  doi          = {10.5281/zenodo.11280343},
  url          = {https://doi.org/10.5281/zenodo.11280343},
}
```

**Dataset DOI:** [10.5281/zenodo.11280343](https://doi.org/10.5281/zenodo.11280343).

## MMDental

**Repository use:** Multimodal CBCT and expert-record availability/provenance in source-only pipelines (for example X57, X68 and X77). The default profiles do not read its medical records.

**Source:** [Figshare dataset record](https://doi.org/10.6084/m9.figshare.28505276.v1).

**Licence and access:** CC BY ([4.0 terms](https://creativecommons.org/licenses/by/4.0/)) in the dataset record; attribution, licence link and disclosure of changes. Public dataset metadata/download links are provided without a stated registration requirement. No data or anatomical derivatives are bundled.

**Exact dataset citation supplied by Figshare:**

> Wang, Chengkai; Zhang, Yifang; Wu, Chengyu; Huang, Xingliang; Wu, Liuxi; Wang, Yitong; et al. (2025). MMDental - A multimodal dataset of tooth CBCT images with expert medical records. figshare. Dataset. https://doi.org/10.6084/m9.figshare.28505276.v1

The dataset licence is **CC BY 4.0**. The [descriptive article](https://doi.org/10.1038/s41597-025-05398-7) is published under **CC BY-NC-ND 4.0**; these are separate objects. Dataset authors in the record differ from the corrected article byline; the record's own citation is preserved.

## Mandibular defects

**Repository use:** X4 source-only mandibular reconstruction geometry. Retain the selected version/member manifest; catalogue geometry counts need not equal the current record's model count.

**Source:** [Figshare dataset record](https://doi.org/10.6084/m9.figshare.28052240.v2).

**Licence and access:** CC BY 4.0 ([4.0 terms](https://creativecommons.org/licenses/by/4.0/)) in the dataset record; attribution, licence link and disclosure of changes. Public dataset metadata/download links are provided without a stated registration requirement. No data or anatomical derivatives are bundled.

**Exact dataset citation supplied by Figshare:**

> Shao, Liangjing (2024). A Mandibular Defect Dataset for Autonomous Reconstruction Planning in Oral and Maxillofacial Surgery. figshare. Dataset. https://doi.org/10.6084/m9.figshare.28052240.v2

## Open-Full-Jaw

**Repository use:** Patient-specific jaw, tooth and periodontal-ligament FE-ready geometry named by the source-only geometry/support implementations (for example X18). No upstream meshes, simulation files or code are bundled as dataset assets.

**Source:** [Authors' repository](https://github.com/diku-dk/Open-Full-Jaw) and [upstream LICENSE](https://github.com/diku-dk/Open-Full-Jaw/blob/main/LICENSE). The repository licence is **CC BY-NC-SA 4.0**: attribution, noncommercial use and ShareAlike for distributed adaptations. Public access is provided; dataset retrieval uses Git LFS. The article's licence is a separate object.

**Citation requested in the upstream README, verbatim:**

```bibtex
@article{gholamalizadeh2022open,
  title     = {Open-Full-Jaw: An open-access dataset and pipeline for finite element models of human jaw},
  author    = {Gholamalizadeh, Torkan and Moshfeghifar, Faezeh and Ferguson, Zachary and Schneider, Teseo and Panozzo, Daniele and Darkner, Sune and Makaremi, Masrour and Chan, Fran{\c{c}}ois and S{\o}ndergaard, Peter Lampel and Erleben, Kenny},
  journal   = {Computer Methods and Programs in Biomedicine},
  volume    = {224},
  pages     = {107009},
  year      = {2022},
  publisher = {Elsevier}
}
```

**DOI:** [10.1016/j.cmpb.2022.107009](https://doi.org/10.1016/j.cmpb.2022.107009).

## Alsheghri / intellident-ai teethPreparationData

**Repository use:** Prepared-incisor mesh and margin/segmentation input for XBREAK_HUNT_5. Its source acquisition script names the authors' `test Data` VTP subset. No VTP file, training/test geometry or learned weight is bundled.

**Source:** [intellident-ai/teethPreparationData](https://github.com/intellident-ai/teethPreparationData) and its [LICENSE](https://github.com/intellident-ai/teethPreparationData/blob/main/LICENSE). **MIT License; Copyright (c) 2024 intellident-ai.** The licence permits use, modification and redistribution subject to retaining its copyright and permission notice. There is no README or mandated scholarly citation in the retrieved repository root. Public access has no stated registration requirement.

**Repository attribution (no upstream-mandated bibliographic format):**

> intellident-ai. (2024). teethPreparationData [Dataset]. https://github.com/intellident-ai/teethPreparationData

The repository is linked by the authors' article [Mesh based segmentation for automated margin line generation on incisors receiving crown treatment](https://arxiv.org/abs/2507.22859), Ammar Alsheghri, Ying Zhang, Farnoosh Ghadiri, Julia Keren, Farida Cheriet and Francois Guibault (2025), DOI [10.48550/arXiv.2507.22859](https://doi.org/10.48550/arXiv.2507.22859). This is supplementary scholarly attribution, not an invented upstream citation requirement. When redistributing any upstream asset, retain the full MIT notice with that asset; Apache-2.0 does not replace it.

## Descriptive dataset articles

These supplementary references explain the external datasets; they are not required citation text invented for the source releases. No article is bundled.

Wang Y, Ye F, Chen Y, Wang C, Wu C, Xu F, Ma Z, Liu Y, Zhang Y, Cao M, Chen X. (2025). A multi-modal dental dataset for semi-supervised deep learning image segmentation. Scientific data 12: 117. [DOI](https://doi.org/10.1038/s41597-024-04306-9). Article: **CC BY-NC-ND 4.0**; the linked Zenodo/Figshare datasets are separately CC BY 4.0.

Wu J, Jiang L, Shao L, Wang W, Xu X, Zhou Y, Wang X, Wang J, Wu J, Chen X, Zhang S. (2025). A Mandibular Defect Dataset for Autonomous Reconstruction Planning in Oral and Maxillofacial Surgery. Scientific data 12: 1763. [DOI](https://doi.org/10.1038/s41597-025-06048-8). Article: **CC BY-NC-ND 4.0**; the linked Zenodo/Figshare datasets are separately CC BY 4.0.

Wang, Chengkai, Zhang, Yifan, Wu, Chengyu, Liu, Jun, Huang, Xingliang, Wu, Liuxi, Wang, Yitong, Feng, Xiang, Lu, Yiting, Wang, Yaqi. (2025). MMDental - A multimodal dataset of tooth CBCT images with expert medical records. Scientific Data 12, 1172. [DOI](https://doi.org/10.1038/s41597-025-05398-7). Article: **CC BY-NC-ND 4.0**; dataset: CC BY 4.0.

Filippon, I.; Tanner, C.; von Jackowski, J.A.; Schulz, G.; Toepper, T.; Müller, B. (2024). Determining Aligner-Induced Tooth Movements in Three Dimensions Using Clinical Data of Two Patients. Oral 4(4), 487–504. [DOI](https://doi.org/10.3390/oral4040039). Article and Zenodo data: CC BY 4.0. The article is linked by the Zenodo record.

## Published numerical table sources

Only numerical facts and source locators are bundled in the small fixtures; articles, source table images and figures are excluded. Retrospective derived contrasts are identified as derived values, not new measurements. Citations below are formatted from primary article metadata (publisher deposits and article XML); unlike the dataset-page blocks above, no fixed citation wording is mandated.

### Cunali et al. (2017)

**Repository use and locator:** X59 and the R4C laboratory contract reproduce regional micro-CT versus silicone-replica group summaries. Table 1; Figure 4, p.470; source methods pp.468–469.

Cunali RS, Saab RC, Correr GM, Cunha LFD, Ornaghi BP, Ritter AV, Gonzaga CC. (2017). Marginal and Internal Adaptation of Zirconia Crowns: A Comparative Study of Assessment Methods. Brazilian dental journal 28(4): 467-473. https://doi.org/10.1590/0103-6440201601531

**Terms:** CC BY 4.0, linked by the [primary SciELO article](https://www.scielo.br/j/bdj/a/SqXx7S5vR7yn3gSCdTkMrrK/?lang=en); attribution, licence link and disclosure of changes. [Primary article](https://doi.org/10.1590/0103-6440201601531).

### Sagheb et al. (2023)

**Repository use and locator:** X63 and R4P reproduce first/tenth tightening preload summaries. Results; Figure 5. Table 1 gives tightening/friction data; it does not contain the endpoint preload means.

Sagheb K, Goergen CI, Doell S, Schmidtmann I, Wentaschek S. (2023). Preload and friction in an implant-abutment-screw complex including a carbon-coated titanium alloy abutment screw: an in vitro study. International journal of implant dentistry 9(1): 8. https://doi.org/10.1186/s40729-023-00473-3

**Terms:** CC BY 4.0; attribution, licence link and disclosure of changes. [Primary article](https://doi.org/10.1186/s40729-023-00473-3).

### Pulpotomy primary study

**Repository use and locator:** X88 reproduces symptom category versus treatment-pathway counts. Table 4.

Patel N, Khan I, Jarad F, Zavattini A, Koller G, Pimentel T, Mahmood K, Mannocci F. (2025). The short-term postoperative pain and impact upon quality of life of pulpotomy and root canal treatment, in teeth with symptoms of irreversible pulpitis: A randomized controlled clinical trial. International endodontic journal 58(1): 55-70. https://doi.org/10.1111/iej.14144

**Terms:** CC BY 4.0. [Primary article](https://doi.org/10.1111/iej.14144).

### Pulpotomy comparator

**Repository use and locator:** Comparator allocation context for X88; the preserved 32/33 denominator ambiguity is not repaired. Chart 1; source allocation description.

Baranwal HC, Mittal N, Yadav J, Rani P, Naveen Kumar PG. (2022). Outcome of partial pulpotomy verses full pulpotomy using biodentine in vital mature permanent molar with clinical symptoms indicative of irreversible pulpitis: A randomized clinical trial. Journal of conservative dentistry : JCD 25(3): 317-323. https://doi.org/10.4103/jcd.jcd_118_22

**Terms:** CC BY-NC-SA 4.0; attribution, noncommercial use and ShareAlike. Only numerical facts and locators are included, no source chart or prose. [Primary article](https://doi.org/10.4103/jcd.jcd_118_22).

### Regional milling trueness

**Repository use and locator:** X89 source context and the milling fixture reproduce regional group RMS means. Table 1.

Son K, Lee JH, Lee KB. (2021). Comparison of Intaglio Surface Trueness of Interim Dental Crowns Fabricated with SLA 3D Printing, DLP 3D Printing, and Milling Technologies. Healthcare (Basel, Switzerland) 9(8): 983. https://doi.org/10.3390/healthcare9080983

**Terms:** CC BY 4.0. Regional milling RMS is not isolated cutter deflection. [Primary article](https://doi.org/10.3390/healthcare9080983).

### Varga et al. (2020)

**Repository use and locator:** X87 uses mandibular implant-position group means and SDs. Published Table 4, MANDIBLE; accepted-manuscript Table 3.

Varga E, Antal M, Major L, Kiscsatári R, Braunitzer G, Piffkó J. (2020). Guidance means accuracy: A randomized clinical trial on freehand versus guided dental implantation. Clinical oral implants research 31(5): 417-430. https://doi.org/10.1111/clr.13578

**Terms:** CC BY-NC 4.0 in the publisher-deposited Crossref licence field (version of record); attribution and noncommercial use. Only numerical facts and locators are included. [Primary article](https://doi.org/10.1111/clr.13578).

### Dynamic-navigation study

**Repository use and locator:** X87 uses dynamic-navigation position group summaries. Results / Comparison of accuracy; Figure 4.

Wu D, Zhou L, Yang J, Zhang B, Lin Y, Chen J, Huang W, Chen Y. (2020). Accuracy of dynamic navigation compared to static surgical guide for dental implant placement. International journal of implant dentistry 6(1): 78. https://doi.org/10.1186/s40729-020-00272-0

**Terms:** CC BY 4.0. [Primary article](https://doi.org/10.1186/s40729-020-00272-0).

### Implant-placement meta-analysis

**Repository use and locator:** Guide fixture retains navigation subgroup context; pooled-mean confidence intervals are not individual-case SDs. Figures 8–10.

Khaohoen A, Powcharoen W, Sornsuwan T, Chaijareenont P, Rungsiyakull C, Rungsiyakull P. (2024). Accuracy of implant placement with computer-aided static, dynamic, and robot-assisted surgery: a systematic review and meta-analysis of clinical trials. BMC oral health 24(1): 359. https://doi.org/10.1186/s12903-024-04033-y

**Terms:** CC BY 4.0; article-data CC0 exception as stated in its permissions unless a credit line says otherwise. [Primary article](https://doi.org/10.1186/s12903-024-04033-y).

### Crown technician-variation study

**Repository use and locator:** CROWN_DIAG retains a same-technician repeated CAD-file RMS reference. Table 1.

Liu CM, Lu TY, Wang CS, Feng SW, Lin YC, Lee SY, Lin WC. (2025). Evaluation of the accuracy, occlusal contact and clinical applications of zirconia crowns using artificial intelligence design versus human design. Journal of dental sciences 20(3): 1665-1672. https://doi.org/10.1016/j.jds.2025.03.016

**Terms:** CC BY 4.0. [Primary article](https://doi.org/10.1016/j.jds.2025.03.016).

### Nagata et al. (2025)

**Repository use and locator:** CROWN_DIAG retains conventional CAD variation at six occlusal cusp positions after manufacture and scan. Table 1 and methods.

Nagata K, Inoue E, Nakashizu T, Seimiya K, Atsumi M, Kimoto K, Kuroda S, Hoshi N. (2025). Verification of the accuracy and design time of crowns designed with artificial intelligence. The journal of advanced prosthodontics 17(1): 1-10. https://doi.org/10.4047/jap.2025.17.1.1

**Terms:** CC BY-NC 4.0; attribution and noncommercial use. Only numerical facts and locators are included, no article text or figures. [Primary article](https://doi.org/10.4047/jap.2025.17.1.1).

### Choi et al. (2019)

**Repository use and locator:** R4F source contract names protocol-specific static and fatigue observations; the printed lab plan requests future same-specimen measurements. Tables 2 and 3; Methods: 30 degrees, 11 mm, R=0.1, 15 Hz, 5 million cycles.

Choi NH, Yoon HI, Kim TH, Park EJ. (2019). Improvement in Fatigue Behavior of Dental Implant Fixtures by Changing Internal Connection Design: An In Vitro Pilot Study. Materials (Basel, Switzerland) 12(19): E3264. https://doi.org/10.3390/ma12193264

**Terms:** CC BY 4.0; no source article is bundled. [Primary article](https://doi.org/10.3390/ma12193264).

## Manufacturer document conditions

The implant profiles retain three numerical depth/datum conditions from manufacturer documents. Manufacturer copyright remains with each source; no document, illustration, drill CAD model or complete table is distributed. These are generic research scenarios, not instructions for a clinical procedure. No scholarly citation wording is mandated by the fixtures.

- **Straumann_BLX_2024_VeloDrill**: 702115, PDF page9 / printed8, 20-Nov-2024. [Source document](https://www.straumann.com/content/dam/media-center/straumann/en/documents/brochure/technical-information/702115-en_low.pdf#page=9). SHA-256 `4532ac919ac2bd557c8d564a8d1b3b7dc2bd3aa971fe695624842c6c9bf3ef38`.
- **Astra_EV_Guided_2017**: 32671186-USX-1712, printed9 and measurement table printed26. [Source document](https://www.dentsplysirona.com/content/dam/dentsply/pim/manufacturer/Implants/Implant_systems/Astra_Tech_Implant_System_EV/Surgical/Accessories/Surgical_accessories_EV___GS/32671186-USX-1712%20Guided%20surgery%20Astra%20Tech%20Implant%20System%20EV_LR.pdf#page=9). SHA-256 `165a88b726ba9c64fda615d45b8067f3834bdfe02b77f8286872c1f612b3698d`.
- **ZimVie_T3_parallel_ACT3p85_datum_aligned**: ZVINST0012-GLBL REV E04/26, Tables printed51 and53; crestal/subcrestal reference printed54. [Source document](https://www.zimvie.eu/content/dam/zimvie-corporate/en/dental/literature/zvinst0012/zvinst0012_surgical_manual_t3-t3_pro_osseotite_final_secured.pdf#page=51). SHA-256 `d4b67b2f2737ac41fc339b7753541912e88d2da5b6ee14a7183a4a7440f7681b`.

Other DOI references in source-only implementations remain attached to their original functions and source locators; they are not additional bundled source articles or a permission to redistribute them.
