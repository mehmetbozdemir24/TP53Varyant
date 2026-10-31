# TP53Varyant Teknik Kılavuz

Bu belge kodu devralacak, sonuçları yeniden üretecek veya yöntemi değerlendirecek kişi içindir. Kullanıcı
gözünden anlatım için [KULLANIM_REHBERI.md](KULLANIM_REHBERI.md) dosyasına bakın. Buradaki tüm eşikler, sayılar
ve dosya adları kodun 31 Ekim 2026 tarihli halinden alınmıştır.

## İçindekiler

1. [Mimari ve veri akışı](#1-mimari-ve-veri-akışı)
2. [Veri kaynakları ve okuyucular](#2-veri-kaynakları-ve-okuyucular)
3. [Birleşik tablo: variants.parquet](#3-birleşik-tablo-variantsparquet)
4. [Kural motoru: ClinGen TP53 VCEP v2.4](#4-kural-motoru-clingen-tp53-vcep-v24)
5. [Franklin-benzeri mod](#5-franklin-benzeri-mod)
6. [Makine öğrenmesi](#6-makine-öğrenmesi)
7. [LLM ile kanıt çıkarımı](#7-llm-ile-kanıt-çıkarımı)
8. [Arayüz (app.py)](#8-arayüz-apppy)
9. [Doğrulama ve test](#9-doğrulama-ve-test)
10. [Bakım ve genişletme](#10-bakım-ve-genişletme)
11. [Bilinen sınırlar](#11-bilinen-sınırlar)
12. [Kaynaklar](#12-kaynaklar)

---

## 1. Mimari ve veri akışı

```
scripts/download.py ──► data/raw/*            (ham dosyalar + MANIFEST.json)
        │
python -m tp53.build ──► data/processed/variants.parquet   (16.136 varyant × 152 kolon)
        │                 data/processed/gene.json, erepo.parquet, nci_*.parquet, civic_*.parquet
        │
        ├── scripts/esm2_llr.py ──► data/processed/esm2_llr.csv      (isteğe bağlı, GPU)
        ├── python -m tp53.ml   ──► data/processed/ml_predictions.parquet, models/*.joblib, models/metrics.json
        └── scripts/validate_*.py ──► data/processed/validation_*.csv
                │
streamlit run app.py ── tp53.acmg (kural motoru, çalışma anında) ── tp53.llm (yalnızca düğmeyle, canlı)
```

`run_pipeline.py` bu adımları sırayla çağırır: `--download` ile indirme, ardından build, VCEP doğrulaması,
`--esm2` ile ESM-2 ve son olarak ML. Tüm hat bir dizüstü bilgisayarda yaklaşık 3-4 dakika sürer (ESM-2 hariç).

Tasarım kararları:

- **Anahtar:** her varyant `NM_000546.6` üzerindeki normalize `c.` notasyonuyla tanımlanır (`tp53.common.norm_c`).
  Kaynakların hepsi bu anahtara bağlanır.
- **Önceden hesaplama:** veri birleştirme, ML skorları ve yüzdelikler önceden hesaplanır. Arayüz çalışırken
  yalnızca kural motoru (milisaniyeler) ve isteğe bağlı canlı sorgular (ClinVar efetch, PubMed, Ollama) çalışır.
- **Karar verici:** sınıf yalnızca kural motorundan çıkar. ML skorları ve LLM çıktısı puana kendiliğinden eklenmez.

| Modül | Satır | Görev |
|---|---|---|
| `tp53/common.py` | 188 | Yollar, transkript sabitleri, HGVS normalizasyonu, sonuç tipi, domainler, Grantham, kodon tablosu |
| `tp53/build.py` | 650+ | Kaynak okuyucular, CDS yeniden kurulumu, birleştirme, `variants.parquet` |
| `tp53/acmg.py` | 780 | `Evidence`, puan → sınıf, VCEP motoru, klinik kriterler, Franklin modu, arama (`find`) |
| `tp53/ml.py` | 228 | Model A ve Model B eğitimi, çapraz doğrulama, kalibrasyon |
| `tp53/llm.py` | 199 | ClinVar VCV XML ayrıştırma, PubMed, Ollama ile yapılandırılmış çıkarım |
| `tp53/help_texts.py` | 280+ | Arayüzdeki açıklama metinleri, araç tanımları, HGVS rehberi |
| `app.py` | 800 | Streamlit arayüzü (9 sekme) |

---

## 2. Veri kaynakları ve okuyucular

Her indirme `data/raw/MANIFEST.json` dosyasına URL, UTC zaman, bayt ve md5 ile yazılır. Tek bir kaynak
yenilemek için: `python scripts/download.py gnomad`. Anahtarlar: `vcep, nci, uniprot, hotspots, erepo,
myvariant, gnomad, clinvar, civic`.

| Kaynak | Uç nokta | Dosya | Okuyucu |
|---|---|---|---|
| ClinGen TP53 VCEP Tablo S1-S3 | Springer ESM: `13073_2025_1536_MOESM{1,2,3}_ESM.xlsx` | `vcep_TableS{1,2,3}.xlsx` | `read_vcep` |
| MyVariant.info (dbNSFP, CADD, dbscSNV) | `myvariant.info/v1/query?q=dbnsfp.genename:TP53&fetch_all=true` (scroll) | `myvariant_TP53_dbnsfp.json` | `read_dbnsfp` |
| gnomAD v4.1 | GraphQL `region(chrom 17, 7661000-7688500, GRCh38)`, `gnomad_r4` | `gnomad_v4_TP53.json` | `read_gnomad` |
| NCI TP53 Database r21 | `tp53.cancer.gov/static/data/{Tumor,Germline,Function}Download_r21.csv` | `nci_*_r21.csv` | `read_nci` |
| ClinVar | E-utilities `esearch` (`TP53[gene] AND single_gene[prop]`) + `esummary` (300'lük parçalar) | `clinvar_TP53_esummary.json` | `read_clinvar` |
| ClinGen eRepo | `erepo.clinicalgenome.org/evrepo/api/classifications?gene=TP53&format=json&matchLimit=5000` | `erepo_TP53.json` | `read_erepo` |
| UniProt P04637 | `rest.uniprot.org/uniprotkb/P04637.json` ve EBI Proteins `variation/P04637` | `uniprot_*.json` | `read_uniprot` |
| cancerhotspots.org v2 | `cancerhotspots.org/api/hotspots/single` (TP53 satırları) | `cancerhotspots_single_TP53.json` | `read_hotspots` |
| CIViC | `civicdb.org/downloads/nightly/nightly-{ClinicalEvidenceSummaries,VariantSummaries}.tsv` | `civic_*_TP53.tsv` | `read_civic` |

### 2.1 Kaynak bazında ayrıntılar ve tuzaklar

**VCEP Tablo S1 (PVS1, kanonik splice).** Tablo, her donör/akseptör bloğunda (bir "WT" satırı ve altındaki
6 varyant) gücü birleştirilmiş hücrede verir. pandas bunu yalnızca ilk satıra okur. Çözüm: "WT" ile başlayan
satırlarda `cumsum` ile blok numarası verilip 15, 16 ve 13. kolonlar blok içinde `ffill/bfill` ile doldurulur.
Çıkan kolonlar: `s1_pvs1`, `s1_rationale`, `s1_rna`, `s1_pred1`.

**VCEP Tablo S2 (PP3/BP4).** Tüm missense için Align-GVGD sınıfı, BayesDel (no-AF) ve panelin kodu
(`s2_code`: `PP3_moderate`, `PP3`, `BP4`, `BP4_moderate`, `No evidence`). Kolonlar: `p3, agvgd,
bayesdel_vcep, s2_code, s2_spliceai, s2_splice_note`. S2'de olmayan missense için `bayesdel_vcep`, dbNSFP
BayesDel ile doldurulur ve kod motorda yeniden hesaplanır (bkz. 4.3).

**VCEP Tablo S3 (PS3/BS3).** Kato, Funk, Giacomelli, Kotler, Kawaguchi sınıfları ve panel kodu (`s3_code`).
Aynı protein değişimi farklı nükleotid değişimleri için ayrı satır olabilir. Birleştirmede `s3_seq_change`
ile c. değişimi eşleşen satır tutulur, eşleşme yoksa ilki.

**dbNSFP (MyVariant).** Kayıtlar hg19 `_id` taşır ve transkript bazında liste döner. `_agg` yüksek-zararlı
araçlarda en büyüğü, düşük-zararlı araçlarda (SIFT, SIFT4G, FATHMM, PROVEAN) en küçüğü alır. MutationTaster'ın
`score` alanı tahminin doğruluk olasılığıdır; `_mt_oriented`, tahmin `D/A` ise skoru, `N/P` ise `1 - skor`
alarak zararlılık yönüne çevirir. dbNSFP satırları c. anahtarına `build_position_map` ile çıkarılan hg38 → c.
eşlemesi üzerinden bağlanır: 3.318 SNV'nin 2.893'ü eşlenir, 425'i (çoğu intron derinliği) eşlenemez.

**gnomAD v4.1.** Gen sorgusu FAF ve UTR/intron varyantlarını vermediği için bölge sorgusu kullanılır
(16.990 varyant). Exome ve genome sayımları toplanır; alt-grup kimlikleri `_` içerenler ve `XX/XY` atlanır.
FAF95, gnomAD'ın yöntemiyle yerel hesaplanır:

```
FAF95(ac, an) = chi2.ppf(0.05, 2·ac) / 2 / an
```

`gnomad_faf95_popmax` yalnızca kurucu olmayan gruplardan, AN ≥ 2000 ve AC ≥ 2 koşuluyla seçilir. Dışlanan
gruplar (`FOUNDER_POPS`): `asj, fin, ami, mid, remaining`. PM2 için ayrıca `gnomad_max_pop_af_multi` (kurucu
olmayan ve AC > 1 grupların en yüksek AF'si) tutulur. Ham grup sayıları `gnomad_pops` kolonunda JSON olarak durur.

**NCI TP53 r21.** Somatik dosya c. bazında toplanır (`nci_tumor_n` = satır sayısı, diğer kolonlar ilk değer);
SpliceAI dört skoru `spliceai_max` olarak birleştirilir. Tümör topografisi ayrı dosyaya yazılır
(`nci_tumor_topography.parquet`). Germline dosyasından aile, birey ve sınıf sayıları çıkar;
`nci_germ_lfs_families` yalnızca `LFS`, `LFL` ve `TP53_Chompret` sınıflı aileleri sayar (FH, noFH, Other
sayılmaz). Fonksiyonel dosya `nci_function.parquet` olarak saklanır ve p1 bazında deney sayısı eklenir.

**ClinVar.** `title` alanı `NM_000546.x(TP53):c....` deseniyle ayrıştırılır; desene uymayan (çok genli) kayıtlar
atlanır. Germline, onkojenite ve somatik klinik etki ayrı tutulur. Yıldız (`CV_RANK`): practice guideline 4,
expert panel 3, multiple submitters no conflicts 2, single submitter veya conflicting 1, diğerleri 0. Aynı c.
için en yüksek yıldızlı kayıt tutulur. Sonuç: 3.890 kayıt.

**ClinGen eRepo.** `matchLimit` verilmezse API yalnızca 25 kayıt döner; 5000 ile 189 TP53 sınıflaması gelir.
Her kayıttan sınıf, kılavuz sürümü, tarih ve `Met` durumundaki kodlar alınır. `code_points`, kodları VCEP
Bayes puanına çevirir (`PS3_Moderate` → 2, `BS1` → -4, `BA1` → -8).

**UniProt.** Çekirdek kayıttan bölge, motif, bağlanma yeri vb. özellikler (`gene.json`) ve doğal varyantların
UniProt notu (`uniprot_desc`, örn. "in LFS; germline mutation ...") alınır. EBI Proteins varyasyon API'si
ClinVar/Ensembl kaynaklı klinik anlam ekler; inceleme durumu `no assertion criteria provided` (0 yıldız) olanlar
atlanır (`uniprot_sig`). Bu filtre olmadan P72R tek bir değerlendirilmemiş gönderim yüzünden "Pathogenic"
görünüyordu.

**cancerhotspots.org.** Kalıntı başına tümör sayısı, q-değeri, doku dağılımı ve değişim sayıları. Varyanta
`hotspot_residue_tumors` ve tam değişim için `hotspot_aa_count` olarak bağlanır. Indel hotspot dosyası
kullanılmaz.

**CIViC.** Gece dökümünden TP53 geçen satırlar tutulur. Filtre, sütun değerinin tam olarak `TP53` olması veya
`TP53 ` ile başlamasıdır (ilk sürümde alt dize eşleşmesi boş dönüyordu).

### 2.2 CDS yeniden kurulumu ve koordinatlar

Kaynaklardan gelen SNV çiftleri (gnomAD hgvsc + pos, NCI hg38, ClinVar SPDI) ile c. konumu → hg38 konumu ve
referans baz eşlemesi kurulur. Kodlayan konumlardaki boşluklar aynı ekson içinde eksi iplik kuralıyla
(c +1 → g -1) doldurulur. 1.182 nt'lik CDS'de hiçbir kaynakta görülmeyen 9 baz kalır
(c.24, 75, 225, 231, 990, 999, 1092, 1098, 1176). Bunlar dört katlı dejenere kodonların 3. bazıdır; protein
değişmediği için yer tutucu "C" konur. Kurulan CDS çevrildiğinde UniProt dizisiyle birebir aynıdır; build günlüğü
bunu her çalıştırmada denetler.

hg19 konumu `hg38 - 96682` (`HG19_HG38_OFFSET`, TP53 bölgesi için sabit) ile türetilir. Franklin bağlantısı
`chr17-{hg19}-{ref}-{alt}`, gnomAD bağlantısı `17-{hg38}-{ref}-{alt}` biçimindedir.

### 2.3 Sonuç tipi ve protein değişimi

`p1` (tek harfli protein değişimi) öncelik sırası: S2 → gnomAD → ClinVar → NCI → CDS üzerinden çeviri.
`consequence(c, p1)` kuralları, sırayla:

1. Intron ofseti varsa (`c.375+1`): |ofset| ≤ 2 `splice_canonical`, ≤ 8 `splice_region`, aksi `intronic`.
   Bu kontrol UTR'den önce yapılır; aksi halde `c.-29+1G>T` yanlışlıkla 5' UTR sayılıyordu.
2. `c.-` → `5_utr`, `c.*` → `3_utr`.
3. Protein: `^M1(\D|$)` → `start_lost` (yalnızca kodon 1; `M160I` gibi değişimler yakalanmaz), `fs` →
   `frameshift`, sonda `*` → `nonsense`, başta `*` → `stop_lost`, `=` → `synonymous`, harf-sayı-harf →
   `missense`, sonra `inframe_*`.

Çerçeve kayması ve sınıflanamayan indel'ler için `mutate_ptc` değişimi CDS'ye uygular ve yeni stop kodonunu
arar: uzunluk farkı 3'ün katıysa `inframe`, stop bulunursa kodon numarası, CDS içinde stop yoksa `beyond`.
Sonuç `ptc_codon` kolonuna yazılır ve PVS1 gücünü belirler.

---

## 3. Birleşik tablo: variants.parquet

16.136 satır, 152 kolon. Anahtar `c`. Kolonlar kaynağa göre gruplanmıştır:

| Grup | Kolonlar |
|---|---|
| Kimlik ve koordinat | `c, p1, p3, consequence, codon, aa_ref, aa_alt, domain, exon, c_pos, c_offset, hg19_pos, hg38_pos, g_ref, g_alt, rsid, franklin_id, hg38_id, ptc_codon` |
| In-silico (dbNSFP) | `revel, alphamissense, bayesdel, varity_r, varity_er, sift, sift4g, polyphen2_hdiv, polyphen2_hvar, primateai, fathmm, fathmm_xf, mutationtaster, mutationassessor, provean, dann, metalr, metasvm, metarnn, clinpred, vest4, mvp, mpc, esm1b, eve, gmvp, deogen2, mutformer, phactboost, mutpred, mutpred_mech, lrt, gerp, phylop100, phastcons100, genocanyon, fitcons, cadd_phred, dbscsnv_ada, dbscsnv_rf, grantham` |
| Splice | `spliceai` (NCI ve S2'nin en büyüğü), `spliceai_max, s2_spliceai, s2_splice_note` |
| Popülasyon | `gnomad_id, gnomad_ac, gnomad_an, gnomad_af, gnomad_hom, gnomad_faf95_popmax, gnomad_faf95_pop, gnomad_max_pop_af_multi, gnomad_pops (JSON), gnomad_filters, gnomad_flags, in_gnomad, gnomad2_exome_af, gnomad2_genome_af, exac_af` |
| VCEP tabloları | `agvgd, bayesdel_vcep, s2_code, s1_pvs1, s1_rationale, s1_rna, s1_pred1, s3_seq_change, kato, funk, giacomelli, kotler, kawaguchi, s3_code` |
| ClinVar | `clinvar_vcv, clinvar_id, clinvar_germline, clinvar_review, clinvar_stars, clinvar_last_eval, clinvar_traits, clinvar_n_scv, clinvar_oncogenicity, clinvar_somatic_tier, clinvar_spdi, clinvar_mc, p_clinvar` |
| ClinGen eRepo | `erepo_class, erepo_version, erepo_date, erepo_met, erepo_not_met, erepo_points, erepo_uuid, erepo_caid, erepo_url` |
| NCI | `nci_tumor_n, tcga_icgc_genie_count, cosmic_id, nci_hotspot, kato_nci, dne_lof_class, dne_class, structure_function_class, agvgd_nci, sift_class_nci, polyphen_class_nci, nci_germ_families, nci_germ_individuals, nci_germ_classes, nci_germ_lfs_families, nci_functional_assays, ...` |
| Hotspot | `hotspot_codon` (VCEP sabit listesi), `hotspot_residue_tumors, hotspot_aa_count` |
| UniProt | `uniprot_var, uniprot_sig, uniprot_desc` |

Sonuç tipi dağılımı: intronic 7.383, missense 2.786, frameshift 1.993, 3' UTR 1.363, synonymous 633,
5' UTR 400, inframe_indel 312, indel_other 254, splice_region 242, inframe_deletion 209, nonsense 197,
splice_canonical 164, other 78, inframe_delins 59, inframe_insertion 44, stop_lost 10, start_lost 9.

ML çıktısı ayrı dosyadadır ve arayüz yüklerken `c` üzerinden eklenir: `ml_predictions.parquet` →
`ml_path_score, ml_path_pct, ml_B{1,2}_{functional,partial,nonfunctional}`.

Kolon tipleri: karışık tipli `object` kolonlar parquet uyumu için metne çevrilir. Sayısal kullanımda
`pd.to_numeric(..., errors="coerce")` gerekir (acmg içindeki `_f` ve `_s` yardımcıları bunu yapar).

---

## 4. Kural motoru: ClinGen TP53 VCEP v2.4

Kaynak: ClinGen TP53 VCEP CSpec v2.4 (GN009) ve Fortuno ve ark. 2025 tablo S1-S3. Kod: `TP53Classifier.vcep`.

### 4.1 Puanlama

`Evidence` her kriteri tutar: kod, karşılandı mı, güç, gerekçe, kaynak, otomatik mi (`auto`), TP53'te
uygulanabilir mi. Puan `STRENGTH_POINTS`: Very Strong 8, Strong 4, Moderate 2, Supporting 1, Stand Alone 8,
No Influence 0; B ile başlayan kodlarda eksi.

`points_to_class`:

| Koşul | Sınıf |
|---|---|
| BA1 karşılandı | Benign |
| toplam ≥ 10 | Pathogenic |
| 6-9 | Likely Pathogenic |
| ≤ -7 | Benign |
| -2 ile -6 | Likely Benign |
| -1, en az 2 benign kriter ve tek patojenik kriter PM2_Supporting | Likely Benign (VCEP istisnası) |
| diğer | Uncertain Significance |

Sonuç iki kez hesaplanır: tüm kanıtlarla (`points`, `classification`) ve yalnızca `auto=True` kriterlerle
(`auto_points`, `auto_classification`).

### 4.2 PVS1

| Durum | Güç |
|---|---|
| nonsense/frameshift, yeni stop < p.351 (Lys351) | Very Strong (NMD beklenir) |
| yeni stop p.351-355 | Strong |
| yeni stop p.356-393 | Moderate |
| çerçeve kayması doğal stop ötesinde sonlanıyor (`beyond`) | Moderate |
| stop konumu hesaplanamadı, başlangıç kodonu < 330 | Very Strong (varsayım, gerekçede yazar) |
| stop_lost + del/dup/ins (uzama, örn. c.1179dup) | Moderate |
| start_lost | Very Strong |
| splice_canonical, Tablo S1'de var | S1 gücü (Very Strong/Strong/Moderate) veya N/A |
| splice_canonical, S1'de yok | Very Strong (varsayılan) |

`pvs1_full` (Very Strong) ise PS3 kullanılmaz; PVS1 Strong ise PS3 Strong → Moderate'e iner.

### 4.3 PP3 / BP4 / BP7

PVS1 karşılanmışsa PP3 kullanılmaz.

Missense: `s2_code` kullanılır. S2'de yoksa Align-GVGD sınıfı ve BayesDel'den aynı kural yeniden kurulur:

| Koşul | Kod |
|---|---|
| GVGD ≥ C65 ve BayesDel ≥ 0.16 | PP3_Moderate |
| C25-C55 ve BayesDel ≥ 0.16 | PP3 |
| BayesDel ≤ -0.008 ve GVGD < C65 | BP4_Moderate |
| -0.008 < BayesDel < 0.16 ve GVGD < C65 | BP4 |
| diğer | kanıt yok |

SpliceAI ≥ 0.2 ise BP4 uygulanmaz ve PP3 (Supporting) splice gerekçesiyle eklenir (`pp3_splice`).

Synonymous, intronic, splice_region, UTR: SpliceAI ≥ 0.2 → PP3; ≤ 0.1 → BP4. BP7 ayrıca şunu ister:
sinonimde ekson çekirdek motifi dışında olmak (eksonun son 3 bazı ve ilk bazı hariç), intronikte ofset
≥ +7 veya ≤ -21. Tek amino asit delesyonları için VCEP'in ayrı BayesDel tablosu bu sürümde yoktur.

### 4.4 PS3 / BS3

Yalnızca missense ve `s3_code` "No evidence" değilse. PP3 SpliceAI'den geliyorsa protein düzeyi testler
kullanılmaz (VCEP uyarısı). `PS3`, `PS3_Moderate`, `PS3_Supporting`, `BS3`, `BS3_Supporting` kodları doğrudan
alınır; PVS1 ile etkileşim 4.2'deki gibidir.

### 4.5 PM1

Missense için:

| Koşul | Güç |
|---|---|
| kodon ∈ {175, 245, 248, 249, 273, 282} (VCEP sabit listesi) | Moderate |
| cancerhotspots.org'da aynı değişim ≥ 10 somatik gözlem | Moderate |
| 2-9 gözlem | Supporting |

İlk sürümde NCI TCGA/ICGC/GENIE sayısına dayalı bir geri dönüş kuralı vardı; uzman kayıtlarıyla karşılaştırmada
20 yanlış PM1'e karşı 3 doğru verdiği için kaldırıldı. NCI sayısı artık yalnızca gerekçede bilgi olarak yazar.

### 4.6 PM2 / BS1 / BA1

| Kriter | Koşul |
|---|---|
| PM2_Supporting | gnomAD'da yok, veya AF < 0.00003 ve `gnomad_max_pop_af_multi` < 0.00004 |
| BS1 (Strong) | 0.0003 ≤ FAF95 popmax < 0.001 |
| BA1 (Stand Alone) | FAF95 popmax ≥ 0.001 |

VCEP, düşük VAF'lı (< 0.35) gnomAD allellerinin klonal hematopoez olabileceğini ve allel dengesi histogramıyla
düzeltme yapılabileceğini belirtir. Bu düzeltme otomatik değildir; gerekçede not olarak yazar.

### 4.7 PS1 / PM5

Karşılaştırma havuzu ClinGen eRepo uzman kayıtlarıdır (ClinVar değil). SpliceAI ≥ 0.2 ise uygulanmaz.

- PS1: aynı amino asit değişimine sahip başka bir c. değişimi uzman panelince Pathogenic → Strong,
  Likely Pathogenic → Moderate.
- PM5: aynı kodonda farklı missense. Karşılaştırılan varyantın sınıfı PM5 çıkarılarak yeniden hesaplanır
  (`_erepo_class_without_pm5`). Grantham uzaklığı karşılaştırılandan küçükse o varyant sayılmaz. ≥ 2 Pathogenic
  → Strong, 1 Pathogenic → Moderate, klinik kanıtlı (PS2/PS4/PP1) Likely Pathogenic → Supporting.
- PS1 + PM5 toplamı 4 puanı geçemez; geçerse PM5 düşürülür.

### 4.8 Klinik kriterler (`_clinical`)

Öncelik: kullanıcı girdisi > eRepo uzman kodu (`use_expert=True` ve kayıt varsa) > karşılanmadı. Uzman
kodları yalnızca `VCEP_CLINICAL_CODES` için alınır: PS2, PS4, PP1, PP4, BS2, BS4, BP2.

| Kriter | Girdi | Kural |
|---|---|---|
| PS2 | de novo vakalar: (güçlü/orta LFS kanseri, anne-baba doğrulandı mı) | güçlü+doğrulanmış 4, güçlü+varsayılan 2, orta+doğrulanmış 2, orta+varsayılan 1 puan; toplam ≥ 8 Very Strong, ≥ 4 Strong, ≥ 2 Moderate, ≥ 1 Supporting |
| PS4 | Chompret, Classic LFS, HER2+ meme < 40 y | 0.5 × Chompret + 1 × Classic + 0.5 × HER2; ölçek PS2 ile aynı. PM2 karşılanmıyorsa uygulanmaz |
| PP1 | meioz ve aile sayısı | ≥ 7 meioz ve > 1 aile Strong, ≥ 5 Moderate, ≥ 3 Supporting. BA1/BS1 varsa uygulanmaz |
| PP4 | VAF %5-25 ve %25-35 gözlem sayısı | %5-25 ≥ 2 ise Moderate, aksi Supporting. BA1/BS1 varsa uygulanmaz |
| BS2 | 60 yaş üstü kansersiz kadın taşıyıcı (tek kaynak) | ≥ 8 Strong, ≥ 4 Moderate, ≥ 2 Supporting |
| BS4 | segregasyon yok | Strong |
| BP2 | trans/cis P/LP | Supporting (CSpec'te güç net değil) |

Güçlü LFS kanserleri: meme < 31 y, koroid pleksus, adrenokortikal < 18 y, rabdomiyosarkom/osteosarkom < 46 y.
Orta: meme 31-49 y, beyin tümörü < 46 y, ACC 18-49 y, akciğer < 46 y, diğer sarkom < 60 y, hipodiploid ALL,
SHH medulloblastom.

### 4.9 TP53'te uygulanmayanlar

PM3, PM4, PM6, PP2, PP5, BP1, BP3, BP5, BP6: `applicable=False` olarak döner ve gerekçede nedeni yazar
(PM6: tüm de novo kanıt PS2'de; PP5/BP6: ClinGen SVI tarafından kaldırıldı).

---

## 5. Franklin-benzeri mod

`TP53Classifier.franklin` Franklin'in genel ACMG davranışını yaklaşık olarak taklit eder. Kurallar Franklin
yardım merkezinden ve TP53 sayfalarında gözlenen çıktılardan (31 Ekim 2026: R175H, R248Q, Y220C, R273H, R213*,
R337H, P72R, V73M, R110H, R290H) çıkarılmıştır.

| Kriter | Kural |
|---|---|
| PVS1 | VCEP motorunun PVS1 sonucu |
| PM2 | max popülasyon AF < 0.0005 → Moderate |
| BS1 | 0.001 ≤ max AF < 0.05 → Strong |
| BA1 | max AF ≥ 0.05 |
| PS1 / PM5 | ClinVar (≥ 1 yıldız, çelişkisiz P/LP) havuzu; PS1 P varsa Strong, aksi Moderate; PM5 aynı kodonda P varsa Moderate, aksi Supporting |
| PP5 | ClinVar P/LP → No Influence (puan yok) |
| BP6 | ClinVar B/LB (çelişkisiz) → **Strong** (Franklin hâlâ uyguluyor) |
| PM1 | iki benign missense arasındaki pencerede P/LP missense yoğunluğu: n ≥ 30 ve yoğunluk ≥ 0.5 Strong, n ≥ 10 ve ≥ 0.3 Moderate, n ≥ 5 ve ≥ 0.2 Supporting; aynı kodonda benign varsa yok |
| PP2 | gende P/LP missense ≥ 8 ve ≥ 3 × B/LB missense |
| PP3 / BP4 | REVEL ≥ 0.932 Strong, ≥ 0.773 Moderate, ≥ 0.644 Supporting (PM1 varsa bir basamak düşer); REVEL < 0.15 BP4; splice için SpliceAI > 0.5 veya dbscSNV Ada > 0.7 |
| BP7 | sinonim ve SpliceAI < 0.1 |
| PS3 / BS3 | önce VCEP Tablo S3 kodu ve gücü; yoksa Kato (fonksiyonsuz PS3 Strong, kısmen PS3 Supporting, fonksiyonel BS3 Strong) |
| BS2 | eRepo uzman kaydındaki BS2 ve gücü |
| PS4 | NCI'de LFS/LFL/Chompret aile sayısı ≥ 5 Moderate, ≥ 2 Supporting; BA1/BS1 varsa uygulanmaz |

Birleştirme Tavtigian puan sistemiyle yapılır (`points_to_class`). Franklin çelişen kanıtları puanla birleştirir:
R110H sayfasında PM2 + PM5 ile BP6 + BS3_Supporting + BS2_Supporting birlikte "Likely Benign" verir; Richards
2015 tablosu bu durumda VUS verirdi. Richards sonucu karşılaştırma için `richards` alanında döner.

Uyarı: BP6 ClinVar kararını doğrudan kanıt saydığı için Franklin modunun ClinVar ile uyumu kısmen döngüseldir.

---

## 6. Makine öğrenmesi

Kod: `tp53/ml.py`, çalıştırma `python -m tp53.ml`. Yalnızca missense (2.786 varyant).

### 6.0 Özetle: ne eğitiyoruz?

İki sınıflandırıcı eğitiyoruz: **Model A** (2 sınıf) ve **Model B** (3 sınıf; B1 ve B2 iki sürüm). REVEL,
AlphaMissense, SIFT gibi araçlar başka grupların modelleridir. Skorları dbNSFP'de her olası değişim için hazır
yayımlanmıştır; biz bu skorları indirip **girdi (özellik)** olarak kullanırız. ESM-2'yi de biz eğitmedik, yalnızca
çalıştırıp skor aldık. Yani Model A, hazır skorların üzerine kurulan TP53'e özel bir üst modeldir (stacking):
hangi araca ne kadar güvenileceğini TP53'te cevabı bilinen varyantlardan öğrenir.

| | Model A | Model B |
|---|---|---|
| Tür | ikili sınıflandırıcı (lojistik regresyon) | 3 sınıflı sınıflandırıcı (gradient boosting) |
| Girdi | 35 sayı: 32 dbNSFP skoru + Align-GVGD (Tablo S2) + 2 ESM-2 skoru | A'nın 35 sayısı + 24 amino asit/konum özelliği (+ B2'de 4 laboratuvar sonucu, Tablo S3) |
| Etiket | 1 = P/LP, 0 = B/LB (uzman paneli, yoksa ClinVar) | Kato testi: 0 fonksiyonel, 1 kısmen, 2 fonksiyonsuz |
| Eğitim satırı | 407 | 2.670 |
| Çıktı | patojenik olma olasılığı | üç sınıfın olasılığı |

Somut satır, R175H (c.524G>A):

```
Model A girdisi (35):
  revel 0.922 · alphamissense 0.986 · bayesdel 0.546 · varity_r 0.951 · varity_er 0.867 · sift 0 · sift4g 0
  polyphen2_hdiv 1 · polyphen2_hvar 0.999 · primateai 0.655 · fathmm -7.55 · fathmm_xf 0.9 · mutationtaster 1
  mutationassessor 2.58 · provean -4.87 · dann 0.999 · metalr 0.996 · metasvm 0.92 · metarnn 0.991
  clinpred 0.993 · vest4 0.938 · mvp 0.997 · mpc 1.1 · esm1b -10.4 · eve 0.502 · gmvp 0.785 · deogen2 0.996
  cadd_phred 25.1 · gerp 5.41 · phylop100 6.16 · phastcons100 1 · lrt 8e-06      (32 sayı: dbNSFP)
  agvgd_num 25                                                                   (VCEP Tablo S2)
  esm2_llr -5.97 · esm2_entropy 0.122                                            (ESM-2, bizim çalıştırmamız)
Model A etiketi: 1 (ClinGen uzman paneli: Pathogenic)        Model A çıktısı: 0.97

Model B ek girdisi (28):
  codon 175 · domain = DNA bağlanma · hyd -4.5 → -3.2 (Δ +1.3) · hacim 173 → 153 (Δ -20.2)
  yük 1 → 0.1 (Δ -0.9) · grantham 29 · to_pro 0 · from_gly 0 · to_gly 0 · cys_change 0
  zinc_site 0 · dna_contact 0                                                    (bizim kodumuz, structural())
  giacomelli_lof 1 · funk_lof 1 · kotler_lof 1 · kawaguchi yok                   (yalnızca B2; Tablo S3)
Model B etiketi: 2 (Kato: Non-functional)     B1 çıktısı: fonksiyonel 0.01 / kısmen 0.32 / fonksiyonsuz 0.67
```

Model A ve B'nin çıktıları ACMG puanına eklenmez ve LLM'e girdi olarak verilmez; yalnızca arayüzde gösterilir.

### 6.1 Özellikler

- **In-silico (35):** `revel, alphamissense, bayesdel, varity_r, varity_er, sift, sift4g, polyphen2_hdiv,
  polyphen2_hvar, primateai, fathmm, fathmm_xf, mutationtaster, mutationassessor, provean, dann, metalr, metasvm,
  metarnn, clinpred, vest4, mvp, mpc, esm1b, eve, gmvp, deogen2, cadd_phred, gerp, phylop100, phastcons100, lrt,
  agvgd_num` + `esm2_llr, esm2_entropy` (`esm2_llr.csv` varsa).
- **Yapı (24, yalnızca Model B):** kodon, 7 domain göstergesi, hidrofobiklik/hacim/yük için yabanıl, mutant ve
  fark (9), Grantham, prolina dönüşüm, glisinden/glisine dönüşüm, sistein değişimi, çinko bağlama kodonu
  (176, 179, 238, 242), DNA temas kodonu (120, 241, 248, 273, 276, 277, 280, 283).
- **Laboratuvar (yalnızca B2):** Giacomelli, Funk, Kotler LOF ikili göstergeleri.

Etiketli kümede kapsama < %50 (A) veya < %30 (B) olan özellikler atılır. Eksikler eğitim içinde medyanla doldurulur.

### 6.2 Model A: patojenite

| | |
|---|---|
| Etiket | 1 = P/LP, 0 = B/LB. Öncelik: eRepo > ClinVar (≥ 1 yıldız, çelişkisiz) |
| Örnek | 407 (248 P/LP, 159 B/LB; 136 eRepo, 271 ClinVar) |
| Aday modeller | `SimpleImputer(median) → StandardScaler → LogisticRegression(C=0.5)` ve `HistGradientBoosting(max_iter 300, lr 0.05, max_leaf_nodes 15, l2 1.0, min_samples_leaf 10)` |
| Seçim | kodon gruplu 5 katlı `StratifiedGroupKFold`, 3 farklı tohum; ortalama ROC-AUC en yüksek olan |
| Kalibrasyon | `CalibratedClassifierCV(method="sigmoid", cv=5)` |
| Sızıntısız skor | etiketli varyantlar için skor, aynı kalibre boru hattının o varyantı görmediği kattaki tahminidir (OOF) |
| Çıktı | `ml_path_score` (0-1), `ml_path_pct` (missense içindeki sıra) |

İzotonik kalibrasyon ilk sürümde kullanılıyordu; skorları basamaklara yığdığı için (933 varyant tam 0.0) sigmoid'e
geçildi. OOF skorları da artık aynı kalibre ölçekte üretilir.

Sonuçlar (kodon gruplu CV):

| Yöntem | ROC-AUC |
|---|---|
| TP53 ML, lojistik (seçilen) | 0.983 (PR-AUC 0.988) |
| TP53 ML, gradient boosting | 0.977 |
| AlphaMissense (tek başına) | 0.990 |
| BayesDel | 0.978 |
| VARITY_R | 0.972 |
| REVEL | 0.970 |
| ESM1b | 0.963 |
| EVE | 0.955 |
| PolyPhen-2 HVAR | 0.950 |
| ESM-2 LLR | 0.925 |
| CADD | 0.913 |
| PrimateAI | 0.893 |
| SIFT | 0.876 |

Tek araç AUC'leri yalnızca skoru olan alt kümede hesaplanır (kapsam ≈ %95). REVEL, VARITY, ClinPred gibi
araçların ClinVar ile eğitildiği unutulmamalıdır (döngüsellik).

### 6.3 Model B: fonksiyon / hipomorf

| | |
|---|---|
| Etiket | Kato 2003 sınıfı (S3): Functional 0, Partially functional 1, Non-functional 2 |
| Örnek | 2.670 (1.400 / 681 / 589) |
| Model | `HistGradientBoostingClassifier`, 3 sınıf; kodon gruplu 5 katlı CV; etiketli varyantlar için OOF olasılık |
| Çıktı | `ml_B{1,2}_functional`, `_partial`, `_nonfunctional` |

| Model | macro-F1 | Dengeli doğruluk | AUC kısmen vs diğer | AUC fonksiyonsuz vs diğer |
|---|---|---|---|---|
| B1 (in-silico + yapı) | 0.655 | 0.651 | 0.730 | 0.905 |
| B2 (+ laboratuvar) | 0.697 | 0.688 | 0.769 | 0.920 |

Karışıklık matrisi (B1; satır gerçek F/K/Fsuz, sütun tahmin): [[1192, 166, 42], [315, 262, 104], [71, 96, 422]].
Kısmen fonksiyonel sınıfın çoğu "fonksiyonel" tahmin edilir; geliştirme için en uygun hedef budur (bkz. 10.4).

### 6.4 ESM-2

`scripts/esm2_llr.py`, `facebook/esm2_t33_650M_UR50D` ile 393 × 19 missense için maskeli marjinal LLR hesaplar:

```
LLR(pos, mut) = log p(mut | dizi, pos maskeli) − log p(wt | dizi, pos maskeli)     (negatif = zararlı)
```

Her pozisyonda tahmin entropisi de yazılır. GPU (fp16) ile RTX 3070 Ti'da birkaç dakika sürer; model ~2.5 GB.
`--model` ile 3B sürüm seçilebilir.

---

## 7. LLM ile kanıt çıkarımı

Kod: `tp53/llm.py`. Model: Ollama üzerinden yerel `gemma4:12b` (ortam değişkeni `TP53_LLM_MODEL`, sunucu
`OLLAMA_HOST`, varsayılan `http://localhost:11434`). İnce ayar yoktur; sıfır atış (zero-shot).

### 7.1 ClinVar gönderimleri

`clinvar_submissions(variation_id)`: `efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=...` ile VCV XML'i
indirir. `is_variationid` parametresi olmadan efetch sıfır gönderim döndürür. Her `ClinicalAssertion` için SCV,
gönderen, sınıf (germline, somatik etki veya onkojenite), değerlendirme tarihi, yorum metinleri (`Comment` ve
`Description` gibi `Attribute` alanları, ilk 6.000 karakter), PubMed kimlikleri, yöntem ve köken alınır.

### 7.2 Metin seçimi ve istem

1. 80 karakterden kısa veya yalnızca "not provided" olan yorumlar atılır.
2. Her yorumda anahtar kelime sayılır (`de novo, segregat, famil, proband, Li-Fraumeni, Chompret, LFS,
   individual, patient, affected, carrier, mosaic, allele fraction, VAF, HER2, functional, transactivation,
   yeast, dominant`) ve en çok geçenler öne alınır.
3. Her yorum 2.500 karakterle kırpılır, toplam 12.000 karaktere kadar eklenir (tipik 8-10 gönderim).

İstem (İngilizce) modele yalnızca açıkça yazanı çıkarmasını, akrabasız probandları saymasını ve her alan için
tanımları verir. Çıktı JSON şemasıyla zorlanır (`format: SCHEMA`):

`de_novo_cases[{cancer, age, parentage_confirmed, quote}], probands_classic_lfs, probands_chompret,
her2_breast_lt40, segregation_meioses, segregation_families, low_vaf_observations,
healthy_female_carriers_over60, functional_studies[{assay, result, quote}], lof_or_gof_notes,
evidence_quotes[], summary_tr`

İstek seçenekleri: `temperature 0, num_ctx 8192, num_predict 1500, think false, keep_alive 10m`, zaman aşımı
900 sn. 8 GB GPU'da bir çıkarım yaklaşık 100 sn sürer. Düşünme modu açıkken ve 16K bağlamda zaman aşımına
düşüyordu.

### 7.3 Kritere çeviri

`to_clinical_suggestion` LLM çıktısını `Clinical` alanlarına çevirir. De novo vakalarda güç kod ile belirlenir:
meme < 31, koroid pleksus, adrenokortikal < 18 (yaş yoksa evet), rabdomiyosarkom/osteosarkom < 46 → güçlü;
diğerleri orta. Arayüzdeki "Öneriyi klinik forma uygula" düğmesi bu değerleri `pending_clinical` olarak
oturuma yazar, `use_expert=False` yapar ve sayfayı yeniden çalıştırır. Hiçbir şey kalıcı kaydedilmez.

### 7.4 Rapor paragrafı

`explain()` LLM'e yalnızca motorun sonucu ve karşılanan kriterlerin listesini verir; 5-7 cümlelik Türkçe yorum
paragrafı ister (sıcaklık 0.2). Puanı etkilemez, ClinVar metnini okumaz.

---

## 8. Arayüz (app.py)

- `load()` (`st.cache_resource`): parquet'ler, ML tahminleri, `gene.json`, ek tablolar, metrikler ve manifest.
  Veri yeniden üretildiyse uygulama yeniden başlatılmalıdır (önbellek süreç ömrü boyunca tutulur).
- Arama `clf.find(q)`: `rs`, `chr17-pos-ref-alt` (hg19 sonra hg38), `c.`, protein (tek veya üç harf, `p.(...)`),
  `R213R` → `R213=`, `M1?` → tüm start_lost; `NM_000546.x(TP53):` önekleri atılır. Birden fazla eşleşmede SNV,
  uzman kaydı ve ClinVar kaydı olan öne alınır; rsID'de gnomAD AF'si yüksek allel önce gelir.
- Klinik form durumu `st.session_state["{c}|{alan}"]` anahtarlarında tutulur; her varyantın formu ayrıdır.
- Açıklama kutuları `explain(key)` ile `help_texts.HELP` sözlüğünden gelir; kenar çubuğundaki
  "Açıklamaları göster" (`SHOW_HELP`) kapatır.
- Tahmin yüzdeliği orta sıra ile hesaplanır: `(P(X < v) + P(X ≤ v)) / 2`, düşük-zararlı araçlarda ters yön.
- Grafikler Altair; bölge görüntüleyici çubuk + domain şeridi + varyant çizgisi + etiket katmanlarıdır.
  Model B grafiğinde `height=alt.Step(30)` kullanılır; sabit yükseklik Streamlit içinde bantları üst üste
  bindiriyordu.

---

## 9. Doğrulama ve test

| Betik | Ne yapar | Çıktı |
|---|---|---|
| `scripts/validate_vcep.py` | 189 eRepo kaydını `use_expert=True` ile yeniden sınıflar; sınıf ve kriter bazında uyum | `validation_vcep_vs_erepo.csv` |
| `scripts/validate_clinvar.py` | ClinVar ≥ 2 yıldız, eRepo dışı, P/LP veya B/LB (879) kayıtları yalnızca otomatik kanıtla ve Franklin moduyla sınıflar | `validation_clinvar2star.csv` |
| `scripts/test_app.py` | Streamlit `AppTest` ile 25 girdi (her sonuç tipi, tüm yazımlar, geçersiz girdiler); tüm sekmeler istisnasız çizilmeli; klinik form aktarımı puanı değiştirmeli | konsol |
| `scripts/test_llm.py` | ClinVar efetch ve Ollama çıkarımını tek varyantta dener | konsol |

Güncel sonuçlar:

- **eRepo (189):** sınıf uyumu %86,8 (164/189). Kriter uyumu: PVS1 %97, PS1 %99, PS3 %99, BS3 %95, PM1 %97,
  PM2 %95, PM5 %93, PP3 %97, BP4 %73, BP7 %99, BA1 %100, BS1 %100. Farkların çoğu uzman kayıtlarının v2.0-2.3
  ile yapılmış olmasından gelir (Funk 2025 testi yoktu; BP4_Moderate etiketi hiç kullanılmamış).
- **ClinVar ≥ 2 yıldız (879), yalnızca otomatik kanıt:**

  | ClinVar | P/LP | VUS | B/LB |
  |---|---|---|---|
  | P/LP (399) | 358 | 41 | 0 |
  | B/LB (480) | 0 | 302 | 178 |

  Ters yönlü hata 0. Franklin modu: P/LP 375/24/0, B/LB 0/0/480 (BP6 nedeniyle döngüsel).
- **Arayüz testi:** 25/25, istisna yok; R175H klinik form aktarımı +12 → +13.

Tekrar üretmek için: `python run_pipeline.py` (veriler indirilmişse) ve ardından
`python scripts/validate_clinvar.py` ile `python scripts/test_app.py`.

---

## 10. Bakım ve genişletme

### 10.1 Veriyi güncelleme

```
python run_pipeline.py --download     # tüm kaynaklar + build + doğrulama + ML
python scripts/validate_clinvar.py
python scripts/test_app.py
```

Ardından uygulamayı yeniden başlatın. Yeni indirmede kontrol edilecekler: build günlüğünde
"ceviri UniProt ile AYNI" satırı, eRepo kayıt sayısının 25 olmaması (matchLimit), gnomAD hatası
(GraphQL şema değişiklikleri), CIViC satır sayısının 0 olmaması.

### 10.2 VCEP yeni sürüm çıkarsa

1. Yeni CSpec'i ve ek tabloları indirin; tablo sütun düzeni değiştiyse `read_vcep` içindeki kolon listelerini
   ve S1 kolon indekslerini (2, 11, 13, 15, 16) güncelleyin.
2. Eşik değişikliklerini `vcep()` içindeki ilgili blokta yapın (PM2/BS1/BA1, PM1 hotspot listesi
   `HOTSPOT_CODONS`, PVS1 kodon sınırları, klinik ölçekler `PS2_POINTS`, `PS2_SCALE`, `PS4_SCALE`).
3. `validate_vcep.py` ile uyumun düşmediğini kontrol edin; eRepo kayıtlarının sürüm alanına göre farkları ayırın.

### 10.3 Yeni bir tahmin aracı eklemek

1. `read_dbnsfp` içine alanı ekleyin (`_agg` yönüne dikkat).
2. Arayüzde göstermek için `app.py` içindeki `TOOLS` listesine (kolon, ad, yön, zararlı eşiği, benign eşiği, not)
   ve `help_texts.TOOL_INFO` sözlüğüne bir cümle ekleyin.
3. ML'de kullanmak için `ml.INSILICO` listesine ekleyin; düşük-zararlıysa `LOWER_IS_DAMAGING` kümesine de.

### 10.4 Önerilen geliştirmeler

- PM2 için gnomAD allel dengesi (VAF < 0.35) düzeltmesi.
- Tek amino asit delesyonları için VCEP BayesDel tablosu ve RNA tabanlı PVS1 (S1 `s1_rna`).
- Model B'ye UniProt/ProtT5 per-residue embedding veya ESM-2 gizli katmanlarından pozisyon özellikleri
  (boyut indirgeme ve kodon gruplu CV şart).
- Model A skorunu ClinGen olabilirlik oranı çerçevesiyle (Pejaver 2022) puana çevirme; şu an puana eklenmez.
- LLM çıktısı için elle doğrulanmış küçük bir test seti ve alan bazında doğruluk ölçümü.
- Hastalık adlarını ClinVar trait kimlikleri (MONDO/OMIM) ile birleştirme; şu an ad normalizasyonu kullanılıyor.

### 10.5 Sık hatalar

| Belirti | Neden | Çözüm |
|---|---|---|
| Arayüz eski sonuçları gösteriyor | `st.cache_resource` | uygulamayı yeniden başlatın |
| eRepo 25 kayıt | `matchLimit` yok | `download.py erepo` |
| LLM sekmesi "Ollama'ya ulaşılamadı" | Ollama kapalı | `ollama serve`, `ollama pull gemma4:12b` |
| LLM zaman aşımı | uzun bağlam veya düşünme modu | `think: false`, `num_ctx 8192` (varsayılan); daha küçük model |
| ClinVar gönderimleri 0 | efetch'te `is_variationid` eksik | `llm.clinvar_submissions` kullanın |
| Windows'ta Türkçe karakterler bozuk | konsol kodlaması | `$env:PYTHONIOENCODING="utf-8"` |
| `git add` uyarısı "LF will be replaced by CRLF" | Windows satır sonu | zararsız; `git config --global core.autocrlf true` |

---

## 11. Bilinen sınırlar

- Uzman kayıtlarının bir kısmı eski CSpec sürümleriyle yapılmıştır; %100 uyum beklenmemelidir.
- PM1 cancerhotspots.org'un anlamlı hotspot listesine dayanır; listede olmayan kalıntılar için sayı eksik kalabilir.
- Franklin modu yaklaşık bir taklittir; Franklin'in PM1 modeli ve ClinVar metin madenciliği kapalı kaynaktır.
- ML etiketlerinde ve bazı in-silico araçlarda ClinVar kaynaklı döngüsellik vardır.
- Çok nükleotidli (delins) missense varyantlarda dbNSFP skoru olmadığı için ML girdileri medyanla doldurulur;
  bu varyantların ML skorları güvenilir değildir.
- Uygulama araştırma ve eğitim amaçlıdır; klinik tanı veya tedavi kararı için kullanılmamalıdır.

---

## 12. Kaynaklar

- Fortuno C, et al. *Genome Medicine* 2025 (TP53 VCEP v2 spesifikasyonları ve Tablo S1-S3). [doi:10.1186/s13073-025-01536-3](https://doi.org/10.1186/s13073-025-01536-3)
- ClinGen TP53 VCEP CSpec v2.4 (GN009): <https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN009>
- Richards S, et al. *Genet Med* 2015;17:405-424.
- Tavtigian SV, et al. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines.
  *Hum Mutat* 2020;41:1734-1737.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification.
  *Am J Hum Genet* 2022;109:2163-2177.
- Kato S, et al. *PNAS* 2003;100:8424-8429.
- Giacomelli AO, et al. *Nat Genet* 2018;50:1381-1387. Kotler E, et al. *Mol Cell* 2018;71:178-190.
- Lin Z, et al. (ESM-2). *Science* 2023;379:1123-1130.
- Franklin yardım merkezi: <https://help.genoox.com/en/collections/2537985-franklin-acmg-automated-classification>
- gnomAD: <https://gnomad.broadinstitute.org> · ClinVar: <https://www.ncbi.nlm.nih.gov/clinvar/> ·
  ClinGen eRepo: <https://erepo.clinicalgenome.org/evrepo/> · NCI TP53: <https://tp53.cancer.gov> ·
  MyVariant.info: <https://myvariant.info> · cancerhotspots.org: <https://www.cancerhotspots.org> ·
  CIViC: <https://civicdb.org>
