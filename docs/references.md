# References and similar studies

A reading list for the team, grouped by topic. Collected from web searches; the
summaries come from search results and abstracts, not from full reads, so check the
papers before citing specifics. Add new items under the right heading (author, year,
link, one line on why it matters).

## Data we use

- Andrews-Goff V, Gales N, Childerhouse S, Laverick S, Polanowski A, Double M. 2023.
  Data from: Australia's east coast humpback whales: satellite tag derived movements.
  Movebank Data Repository, [doi:10.5441/001/1.294](https://datarepository.movebank.org/handle/10255/move.1660) (CC0).
  Also at the [Australian Antarctic Data Centre](https://data.aad.gov.au/metadata/records/AAS_4600_EastAus_humpback_satellite_tracking_data) (doi:10.26179/mnwc-x310).
  Used in `data-raw/tags/` and `data/humpback_satellite_tags.xlsx`.
- [2008/09 southbound tracking from Eden, NSW (data.gov.au)](https://data.gov.au/data/dataset/3ac9c6d0-b0b1-4935-9d66-5b37a6e31d90) and [AAS 2941](https://data.aad.gov.au/metadata/records/fulldisplay/AAS_2941_humpback_whale_satellite_tracks_Eden_2008). Not yet used.

## East Australian humpbacks: movement and migration

- Andrews-Goff et al. 2023. [Satellite tag-derived movements on breeding grounds, feeding grounds and along the migration](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10729012/). Paper describing the tag dataset.
- Andrews-Goff et al. 2018. [Humpback whale migrations to Antarctic summer foraging grounds through the southwest Pacific Ocean](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6098068/). Scientific Reports. 30 tagged whales.
- UQ 2025. [Decades of surveys show whale migration shift](https://news.uq.edu.au/2025-07-decades-surveys-show-whale-migration-shift). Southbound peak moved from early October (2003) to mid-September (2024); linked by the authors to Southern Ocean conditions (sea ice, krill). News item only: find the underlying paper from that page.
- [Coastal fronts utilised by migrating humpback whales (SCU)](https://researchportal.scu.edu.au/esploro/outputs/journalArticle/Coastal-Fronts-Utilized-by-Migrating-Humpback/991013056511602368). Whales favour cooler water and strong temperature fronts along the East Australian Current.
- [Factors impacting migration paths of humpbacks in Queensland (AGU abstract)](https://agu.confex.com/agu/os16/preliminaryview.cgi/Paper90848.html).

## Population recovery (confounds sighting trends)

- Noad et al. 2011. [2004 survey: about 7,090 whales, 10.6% a year growth](https://journal.iwc.int/index.php/jcrm/article/download/318/83). J. Cetacean Res. Manage.
- Paton et al. 2011. [Capture-recapture estimate](https://journal.iwc.int/index.php/jcrm/article/download/317/81).
- UQ 2019. [Reaching peak whale abundance](https://uq.edu.au/news/node/123479). About 10.9% a year growth and about 25,000 whales in 2015; possible slowdown around the mid-2020s.
- [GBRMPA Outlook Report, humpback whales](https://outlookreport.gbrmpa.gov.au/threats-responses-and-risks/8-resilience/83-case-studies-recovery-and-decline-ecosystem/837-humpback-whales). States a possible ~40,000 whales; not a peer-reviewed estimate.

## Habitat suitability and environment

- [Frontiers in Marine Science review, 2021](https://frontiersin.org/articles/10.3389/fmars.2021.720774/full). Feeding-ground humpbacks prefer upwelling, high chlorophyll and fronts.
- [Humpback and minke whale habitat suitability, Atlantic sector of the Southern Ocean (MaxEnt)](https://epic.awi.de/id/eprint/30956/).
- Dalla Rosa et al. 2012. [Modelling BC humpbacks](https://mmru.ubc.ca/wp-content/pdfs/Dalla%20Rosa%20etal%202012%20-%20modeling%20BC%20humpbacks.pdf). Chlorophyll effects were not consistent across models.
- [Oceania breeding-ground habitat projections (warming above ~28 C)](https://repository.up.ac.za/bitstreams/90dc270b-94d2-41fc-bdf8-89316bcb11f4/download). Confirm this is the right paper before citing.

## Presence-only data and observer-effort bias

- [Large whales in the eastern tropical Pacific, MaxEnt bias (Frontiers in Marine Science, 2018)](https://www.frontiersin.org/articles/10.3389/fmars.2018.00419). Predictions can be distorted unless sampling bias is handled.
- Warton et al. 2013. [Modelling observer bias, PLOS ONE](https://ideas.repec.org/a/plo/pone00/0079168.html). Tested on plants, not whales.
- [Combining data types in one spatio-temporal framework (preprint)](https://www.stat.ubc.ca/node/10009). Applied to killer whales; provisional.

## Movement and tag-track methods

- [Meta-analysis of animal movement using state-space models](https://www.academia.edu/95292099/Meta_Analysis_of_Animal_Movement_Using_State_Space_Models). The kind of model behind the smoothed tag positions.
- [Real-time forecasting of animal tracks (Ecology, 2022)](https://openpolar.no/Record/crwiley:10.1002%2Fecy.3718). Predicting future locations.

## Other tagged populations

- Garrigue et al. 2010. [Satellite-tagged humpbacks, New Caledonia](https://www.marinemammals.gov.au/site/assets/files/1071/garrigue-et-al-2010-j-mammal-satellite-tagged-humpbacks-new-caledonia.pdf). Data not located.
- Hauser et al. 2009. [Satellite-monitored humpbacks, Cook Islands](https://whaleresearch.org/wp-content/uploads/2020/12/Movements-of-satellite-monitored-humpbacks-from-the-Cook-Islands.pdf).

## Gaps (nothing found yet)

- Humpback sightings or migration against ENSO / SOI.
- Colwell predictability applied to whales. The statistical frameworks PDF in the repo root covers the method.
