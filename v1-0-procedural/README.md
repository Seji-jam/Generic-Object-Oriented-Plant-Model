# Procedural Process-Based Crop Model

This repository contains the procedural process-based crop model developed by Sajad Jamshidi in Dr. Diane Wang's lab. This effort is part of an ongoing project aimed at developing a base generic plant model adaptable for exploring Genotype by Environment (GxE) scenarios. The model presented here builds upon and significantly improves the previous prototype plant model based on a comprehensive literature review and the evaluation of other crop models conducted during Fall 2023.

## Key Features
- **Procedural Design**: Initially designed as a procedural model to lay the groundwork for future development into an object-oriented framework.
- **Foundation for GxE Exploration**: Serves as a foundational model for studying the interaction between genetic variation and environmental factors.
- **Improvement and Expansion**: Extends the capabilities of the previous version by incorporating findings and methodologies from recent literature and crop model evaluations.

## Model Variables
### Weather and Nitrogen States

| Parameter | Description |
| --- | --- |
| WetIdx | Index indicating the current wetness state |
| NitIdx | Index for nitrogen status |
| WatInp | Water input value |
| NitApp | Nitrogen applied |
| NitFix | Nitrogen fixation rate |

### Plant Types and Growth Factors

| Parameter | Description |
| --- | --- |
| LegTyp | Type of legume (binary, 0 or 1) |
| C3C4Ty | Plant type: C3 or C4 |
| GrowP | Growth parameter or factor |
| SlopeV | Slope value for growth |
| FallR | Fall rate |

### Crop Fraction Values and Nitrogen Concentration

| Parameter | Description |
| --- | --- |
| GrEff | Growth efficiency |
| CarFV | Carbon fraction in vegetative parts |
| YielGV | Yield growth value |
| FatFr | Fat fraction |
| LignFr | Lignin fraction |
| OAcidF | Organic acid fraction |
| MinNF | Minimum nitrogen fraction |
| LeafNC | Nitrogen concentration in living leaves |

### Temperature and Growth Parameters

| Parameter | Description |
| --- | --- |
| BaseTD | Base temperature for development |
| OptTD | Optimal temperature for development |
| CritTD | Critical temperature for development |
| SenTD | Sensitivity to temperature |
| StartP | Start period |
| EndP | End period |
| InitSP | Initial set point |
| LeafWd | Leaf width |
| RootDm | Root depth maximum |
| CropHt | Crop height |
| MaxEH | Maximum expected height |
| DrySI | Dry site index |
| MaxES | Maximum expected size |
| CarbonF | Carbon factor |
| NUpTk | Nitrogen uptake |
| LeafSA | Leaf surface area OR specific area UNKNOWN CHECKING |
| MinSLN | Minimum specific leaf nitrogen |
| MinRCN | Minimum root carbon nitrogen |
| StemNC | Stem nitrogen concentration |
| RootB | Root biomass |
| MaxAJ | Maximal AJ parameter |
| ValN | Nitrogen value parameter |
| ValJ | J value parameter |
| WatUse | Water usage |

### Seed Weight and Nitrogen Content

| Parameter | Description |
| --- | --- |
| SeedWt | Seed weight |
| SeedNC | Seed nitrogen content |
| BulkD | Bulk density |
| MaxHt | Maximum height |
| DevDurV | Development duration vegetative |
| DevDurR | Development duration reproductive |
| SenSc | Senescence scale |

### Soil Moisture and Carbon Parameters

| Parameter | Description |
| --- | --- |
| PanLS | Pan landscape scale |
| SoilTy | Soil type |
| MinWC | Minimum water content |
| FieldC | Field capacity |
| MaxWC | Maximum water content |
| DropM | Drop model |
| DecPM | Decomposition rate for particulate matter |
| RootPM | Root particulate matter |
| BioOM | Biomass organic matter |
| HumusR | Humus rate |
| TotalOC | Total organic carbon |
| BioCH | Biomass carbon humus |
| FrBioC | Fraction of biomass carbon |
| NRatio | Nitrogen ratio |
| AirR | Air ratio |
| SoilSS | Soil seed survival |
| SeedD | Seed depth |
| TillCT | Tillage carbon turnover |
| TillCP | Tillage carbon partition |
| MulF | Mulch factor |

### CO2 and Respiration Rates

| Parameter | Description |
| --- | --- |
| CO2Atm | CO2 atmospheric concentration |
| ResRF | Respiration root factor |
| ResVF | Respiration vegetative factor |
| TempF | Temperature factor |
| CrushF | Crushing factor |
| NitRSh | Nitrogen respiration shift |
| PrePN | Precipitation nitrogen parameter |
| CarbonB | Carbon biomass |
| CarbonX | Carbon excess |
| TempM | Temperature multiplier |

## Development Notes
This version represents an experimental phase in the model's evolution, focusing on developing a procedural framework as a precursor to a more robust object-oriented model. The next phase of development will involve a transition to an object-oriented design, which will undergo extensive testing and refinement.

## Usage and Contribution
This model is released under the MIT license, allowing for free use, modification, and distribution. Contributors and users are encouraged to participate in the ongoing development and enhancement of the model.

## Disclaimer
Please note that this version has not been extensively tested. It represents a trial development phase aimed at establishing a procedural modeling approach. Feedback and contributions are highly welcomed to aid in the transition to an object-oriented model and its subsequent testing and validation.

## Acknowledgments
Special thanks to Dr. Diane Wang's lab for supporting this project and to all contributors who have played a role in the development and refinement of this crop model.

## License
This project is available under the MIT License. For more details, see the LICENSE file in the repository.
