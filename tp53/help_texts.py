"""Arayuzdeki her bolum icin sade Turkce aciklamalar (konuyu bilmeyen kullanici icin).

app.py bu metinleri bolumlerin basinda "Bu bolum ne anlatiyor?" kutusu olarak gosterir.
Kenar cubugundaki "Aciklamalari goster" anahtari ile kapatilabilir.
{v} -> o an acik olan varyantin adi (ornek: P72R).
"""

HELP = {
    # ------------------------------------------------------------------ ACMG
    "acmg_intro": (
        "Bu sekme, **{v}** varyantının hastalık yapıp yapmadığına dair asıl kararı verir. Karar, uluslararası ACMG "
        "kurallarının TP53 uzman panelince (ClinGen VCEP) TP53'e uyarlanmış haliyle verilir.\n\n"
        "**Nasıl çalışır?** Her kanıt (ör. \"toplumda çok nadir\", \"laboratuvarda proteini bozuyor\") bir **kriter** "
        "sayılır ve puan alır: Destekleyici ±1, Orta ±2, Güçlü ±4, Çok güçlü ±8. Patojenik yöndeki kriterler (P ile "
        "başlar) artı, benign yöndekiler (B ile başlar) eksi puandır.\n\n"
        "**Sonuç:** toplam ≥10 Patojenik · 6-9 Olası patojenik · -1..5 VUS (önemi belirsiz) · -2..-6 Olası benign · "
        "≤-7 Benign. Sağdaki küçük kartlar aynı varyantın hasta bilgisi olmadan, Franklin kurallarıyla ve uzman "
        "panelince aldığı sonuçları yan yana gösterir."
    ),
    "acmg_form": (
        "Bazı kriterler varyantın kendisinden değil, **hastalardan ve ailelerden** gelir: varyant kaç Li-Fraumeni "
        "ailesinde görüldü, ailede kanserle birlikte mi geçiyor, anne-babada yokken çocukta mı ortaya çıktı (de novo). "
        "Bunlar veri tabanından otomatik hesaplanamaz. Elinde bu bilgi varsa buraya gir; puan anında yeniden "
        "hesaplanır. Varsayılan olarak, varyant uzman panelinde kayıtlıysa uzmanların klinik kodları kullanılır."
    ),
    "acmg_cards": (
        "Kararın nasıl çıktığını iki parça gösterir.\n\n"
        "**1. Puan hesabı tablosu:** varyant hakkında bulunan her kanıt bir satırdır. \"Ne diyor?\" kanıtı düz "
        "Türkçe anlatır; Puan sütunu o kanıtın getirdiği puandır (🔴 patojenik yönde artı, 🟢 benign yönde eksi). "
        "Altındaki satır toplamı ve sınıfı gösterir.\n\n"
        "**2. Kanıt grupları:** ACMG kanıtları 12 soruya ayırır (ör. \"Laboratuvar testleri\": protein bozuk mu "
        "çıktı?). Her grubun altında soru yazar; renkli kutu o soruya verilen cevaptır: kutunun adı kriter kodu "
        "(PS3, PM1...), altında gücü ve puanı. Kod okuması: **P** = patojenik, **B** = benign yönde; ikinci harf gücü "
        "söyler: **VS** çok güçlü (8), **S** güçlü (4), **M** orta (2), **P** destekleyici (1). \"_Moderate\" gibi ek, "
        "gücün kuralın varsayılanından farklı olduğunu söyler. Varsayılan olarak sadece kanıt bulunan gruplar "
        "gösterilir.\n\n"
        "\"Kanıt görünümü\" ile aynı varyantı TP53 uzman kurallarıyla veya Franklin'in genel kurallarıyla görebilirsin."
    ),
    # ------------------------------------------------------------------ Varyant degerlendirme
    "va_metrics": (
        "Varyantın kimlik kartı. **Sonuç tipi**: DNA değişikliğinin proteine etkisi (missense = bir amino asit "
        "değişir, nonsense = protein erken kesilir, frameshift = okuma çerçevesi kayar, splice = ekson birleşmesi "
        "bozulur, synonymous = amino asit değişmez). **ClinVar**: laboratuvarların ortak kararı. **gnomAD AC/AN**: "
        "sağlıklı ~800 bin kişide bu alleli taşıyan kromozom sayısı / bakılan kromozom sayısı. **NCI somatik tümör**: "
        "kaç tümör örneğinde görüldüğü. **TP53 ML skoru**: bizim modelimizin 0 (zararsız) – 1 (zararlı) tahmini."
    ),
    "va_region": (
        "TP53 proteininin 393 amino asidi soldan sağa dizilir; her çubuk o pozisyonda sınıflanmış varyant sayısıdır "
        "(kırmızı patojenik, gri belirsiz, yeşil benign). Alttaki renkli şerit proteinin bölgeleridir; turuncu uzun "
        "parça DNA'ya bağlanan çekirdek bölgedir. **Mor kesik çizgi {v}'nin yeridir.** Çizgi kırmızı çubukların "
        "yoğun olduğu bir yerdeyse, varyant zararlı varyantların kümelendiği bir bölgededir."
    ),
    "va_clinical": (
        "Üç bağımsız kaynağın bu varyant hakkında söyledikleri. **ClinGen eRepo**: TP53 uzman panelinin resmi kararı "
        "(en güvenilir). **ClinVar**: laboratuvarların gönderdiği kararlar; parantez içindeki inceleme düzeyi ne kadar "
        "güvenilir olduğunu söyler (\"reviewed by expert panel\" en yüksek). **UniProt**: proteinin referans kaydındaki "
        "not ve NCI'nin Li-Fraumeni ailesi kayıtları."
    ),
    "va_pred": (
        "Her satır bir bilgisayar programının (in-silico araç) \"bu amino asit değişimi proteine zarar verir mi?\" "
        "tahminidir. Araçlar hastaya bakmaz; amino asidin canlılar arasında ne kadar korunduğuna, kimyasal "
        "benzerliğe ve protein yapısına bakar.\n\n"
        "**Çubuğun boyu = sıralama.** TP53'teki 2.786 missense varyant o aracın skoruna göre sıralanır. "
        "Çubuk 90'daysa {v}, TP53 varyantlarının %90'ından daha zararlı bulunmuş demektir; 10'daysa en zararsız "
        "%10'un içindedir.\n\n"
        "**Çubuğun rengi = aracın kendi kararı** (yayımlanmış eşiğine göre): kırmızı zararlı, gri kararsız, yeşil "
        "zararsız. Boy ve renk bazen çelişebilir; o yüzden tek bir araca değil, çoğunluğa bakılır.\n\n"
        "**Önemli:** bu araçların hiçbiri tek başına karar vermez. TP53 uzman paneli sadece BayesDel + Align-GVGD "
        "sonucunu kabul eder (grafiğin üstündeki gri satır); ACMG puanına giren yalnızca odur. Fareyi bir çubuğun "
        "üzerine getirince aracın ne ölçtüğü yazar; tablonun \"Ne ölçer?\" sütununda da var."
    ),
    "va_pop": (
        "Varyantın sağlıklı insanlarda ne sıklıkta görüldüğü (gnomAD v4.1, ~800 bin kişi). Mantık basit: Li-Fraumeni "
        "nadir bir hastalıktır; bir varyant toplumda çok sık görülüyorsa bu hastalığın sebebi olamaz.\n\n"
        "**AC** = varyantı taşıyan kromozom sayısı, **AN** = bakılan kromozom sayısı, **AF** = AC/AN. "
        "7.48e-1 = 0,748 (%74,8); 3.0e-5 = 0,00003 (100 bin kromozomda 3). Yıldızlı (*) gruplar kurucu veya karışık "
        "gruplardır, kurallarda hesaba katılmaz. Sağdaki kartta ACMG eşikleri var: binde 1'den sık (BA1) ise tek başına "
        "benign; 100 bin kromozomda 3'ten nadir (PM2) ise patojenik yönde +1."
    ),
    # ------------------------------------------------------------------ Fonksiyon & ML
    "fn_tests": (
        "Laboratuvarda gerçekten ölçülmüş sonuçlar. Varyantlı p53 proteini maya veya insan hücrelerinde üretilir ve "
        "işini (diğer genleri açıp hücre bölünmesini durdurmak) yapıp yapmadığına bakılır. **Fonksiyonel** = normal "
        "çalışıyor (benign yönde kanıt), **fonksiyonsuz / LOF** = bozuk (patojenik yönde kanıt), **kısmen** = arada. "
        "\"—\" o testte bu varyant için ölçüm olmadığı anlamına gelir. Alttaki \"VCEP ön-kodu\", uzman panelinin bu "
        "testleri birleştirdiği karardır ve ACMG puanına giren odur (PS3 +4 / BS3 -4)."
    ),
    "ml_a": (
        "**Model A** (sınıflandırıcı): REVEL, AlphaMissense, SIFT gibi başka grupların hazır modellerinin bu varyant "
        "için ürettiği 35 skoru girdi olarak alır ve tek bir olasılığa çevirir (Franklin'deki \"Aggregated Prediction\" "
        "fikri). Örnek: R175H için girdi REVEL 0.922, AlphaMissense 0.986, SIFT 0 ... (35 sayı), çıktı 0.97. Doğru cevabı bilinen 407 TP53 varyantıyla (248 patojenik, 159 benign) eğitildi. Çıktı 0-1 arası "
        "\"patojenik olma olasılığı\"dır: 0.9 diyen varyantların yaklaşık %90'ı gerçekten patojeniktir. Eğitim setindeki "
        "varyantlar için, modelin o varyantı hiç görmediği turdaki skor gösterilir. **ACMG puanına eklenmez.**"
    ),
    "ml_b": (
        "**Model B** (3 sınıflı sınıflandırıcı): Model A'nın 35 skoruna amino asit değişiminin özelliklerini (konum, "
        "hacim, yük) ekler ve proteinin ne kadar çalıştığını tahmin eder (Kato maya testi sınıfları): fonksiyonel, kısmen "
        "fonksiyonel (hipomorf), fonksiyonsuz. Üç olasılığın toplamı 1'dir. Kısmen fonksiyonel varyantlar daha hafif "
        "etkili, kanser riskini daha az artıran varyantlardır; laboratuvarda test edilmemiş varyantlar için ipucu verir. "
        "B1 sadece dizi ve tahmin araçlarını, B2 ek olarak diğer laboratuvar testlerini kullanır."
    ),
    "ml_perf": (
        "**ROC-AUC**: rastgele bir patojenik ve bir benign varyant seçilince modelin patojeniğe daha yüksek skor verme "
        "olasılığı (1 kusursuz, 0.5 yazı-tura). **Kapsam**: araç skorunun olduğu varyant oranı. Karışıklık tablosunda "
        "satırlar gerçek sınıf, sütunlar tahmindir; köşegen doğru tahminlerdir."
    ),
    # ------------------------------------------------------------------ Diger sekmeler
    "pubs": (
        "Bu varyantın geçtiği yayınlar. Üstteki tablo NCI TP53 veri tabanındaki Li-Fraumeni ailesi kayıtlarıdır "
        "(her satır bir kişi, kaynak makalesiyle). \"PubMed'de ara\" düğmesi varyantın adıyla PubMed'de canlı arama "
        "yapar (internet gerekir)."
    ),
    "gene_dist": (
        "Bu sekme tek bir varyanta değil, **TP53 geninin tamamına** bakar ve her varyantta aynıdır. Grafik, ClinVar'da "
        "kayıtlı TP53 varyantlarını değişiklik tipine göre gruplar ve her grupta kaçının patojenik / belirsiz / benign "
        "olduğunu gösterir. Örnek okuma: frameshift ve nonsense neredeyse hep kırmızı (proteini bozar), synonymous "
        "neredeyse hep yeşil, missense ise çoğunlukla gri (belirsiz) — tahmin araçları en çok missense için gerekir."
    ),
    "gene_hotspot": (
        "Kanser hücrelerinde TP53'ün en sık mutasyona uğrayan pozisyonları (cancerhotspots.org, ~25 bin tümör). "
        "Çubuk = o pozisyonda mutasyon bulunan tümör sayısı. R273, R248, R175 tüm kanserlerdeki en ünlü TP53 "
        "mutasyonlarıdır. Bir varyant bu pozisyonlardaysa ACMG'de PM1 (+2) kriteri karşılanabilir."
    ),
    "gene_uniprot": (
        "UniProt, proteinlerin uluslararası referans kaydıdır (P04637 = insan p53). Tablo proteinin hangi parçasının ne "
        "iş yaptığını söyler: **DNA binding** DNA'ya tutunan bölge, **Region** başka proteinlerle etkileşen parçalar, "
        "**Motif** kısa tanınmış diziler (ör. çekirdeğe giriş sinyali), **Binding site** çinko bağlayan amino asitler. "
        "start/end = başlangıç ve bitiş amino asidi."
    ),
    "dis_intro": (
        "Hastalık adları ClinVar ve NCI veri tabanlarından gelir; elle yazılmaz. Üstteki etiketler, laboratuvarların "
        "{v}'yi **hangi hastalık açısından değerlendirdiğidir** — varyantın o hastalığa yol açtığı anlamına gelmez "
        "(sonuç ACMG sekmesindedir). Alttaki tablo tüm TP53 varyantlarında en sık geçen hastalıklardır. TP53'ün ana "
        "hastalığı Li-Fraumeni sendromudur: tek kopya yeterli, erken yaşta çok çeşitli kanserlere yatkınlık."
    ),
    "som_intro": (
        "Buraya kadar her şey **germline** (kişinin doğuştan taşıdığı) varyant içindi. Bu sekme **somatik** soruya "
        "cevap verir: bu değişiklik bir tümörün içinde bulunursa ne anlama gelir? Tümör veri tabanlarında kaç kez "
        "görüldüğü, kanser hotspot'u olup olmadığı ve tedavi / hastalık seyri ile ilgili literatür kanıtları. "
        "**ACMG puanına girmez ve tedavi önerisi değildir.**"
    ),
    "som_civic": (
        "CIViC, kanser varyantlarının klinik önemini makalelerden özetleyen açık bir veri tabanıdır; her satır bir "
        "makaleden bir iddiadır. **evidence_level**: A klinikte kabul görmüş · B hasta verisi · C vaka raporu · D "
        "laboratuvar · E dolaylı. **evidence_type**: Prognostic (hastalığın seyri), Predictive (tedaviye yanıt), "
        "Diagnostic, Predisposing. **significance**: iddia edilen etki (Poor Outcome = kötü seyir, Sensitivity = "
        "tedaviye yanıt, Resistance = direnç). **evidence_direction**: Supports = çalışma iddiayı destekliyor, Does Not "
        "Support = etki bulunamadı. **kapsam**: Varyant = tam bu değişiklik, Gen = genel olarak \"TP53 mutasyonu\"."
    ),
    "llm_intro": (
        "Bu sekme, ACMG'deki klinik kanıt formunu doldurmana yardım eder. Üç adım vardır ve **yalnızca 2. adımda yapay "
        "zekâ çalışır**:\n\n"
        "1. **ClinVar gönderimlerini getir**: NCBI'dan bu varyant için laboratuvarların yazdığı yorumları indirir "
        "(yapay zekâ yok, birkaç saniye).\n"
        "2. **LLM ile klinik kanıt çıkar**: en bilgili yorumları bilgisayarındaki yerel dil modeline (gemma) verir; "
        "model metinde açıkça yazan aile, hasta ve de novo sayılarını alıntılarıyla çıkarır (~100 sn). Sınıflandırma "
        "kararı vermez.\n"
        "3. **Öneriyi klinik forma uygula**: bulunan sayıları ACMG sekmesindeki forma kopyalar ve puan yeniden "
        "hesaplanır. Basmadan önce alıntıları kontrol et; hiçbir şey kalıcı kaydedilmez."
    ),
    "data_intro": (
        "Uygulamanın kullandığı tüm veri setleri, nereden ve ne zaman indirildikleri. Altta, uygulamanın kararlarının "
        "TP53 uzman panelinin elle verdiği 189 kararla ne kadar örtüştüğü (doğrulama) yer alır."
    ),
}

# Tahmin araclari: ne olcer (sade Turkce, bir cumle)
TOOL_INFO = {
    "ml_path_score": "Bizim modelimiz: aşağıdaki araçların skorlarını TP53'e özel olarak birleştirir (0-1).",
    "revel": "13 aracı birleştiren meta skor; ClinGen'in PP3/BP4 için önerdiği başlıca araç.",
    "alphamissense": "Google DeepMind'ın protein yapısı ve evrimden öğrenen yapay zekâ modeli.",
    "bayesdel_vcep": "Çok aracı birleştiren skor; TP53 uzman panelinin kabul ettiği iki araçtan biri.",
    "varity_r": "Nadir varyantlara odaklanan, ClinVar ile eğitilmiş makine öğrenmesi modeli.",
    "esm1b": "Milyonlarca protein dizisinden öğrenen dil modeli; değişimi ne kadar 'beklenmedik' buluyor (düşük = zararlı).",
    "eve": "Canlılar arasındaki dizi değişkenliğinden öğrenen derin üretici model.",
    "primateai": "Primatlarda sorunsuz görülen değişimlerden öğrenen derin öğrenme modeli.",
    "cadd_phred": "Genom genelinde 60+ özelliği birleştiren skor; 20 = en zararlı %1, 30 = en zararlı %0,1.",
    "polyphen2_hvar": "Amino asit korunması ve 3B yapıya bakan klasik araç (1'e yakın = zararlı).",
    "sift": "Amino asidin türler arasında korunup korunmadığına bakar (0,05 altı = zararlı).",
    "fathmm": "Korunmuşluk tabanlı; kanser genlerinde neredeyse her değişime 'zararlı' der.",
    "mutationassessor": "Protein ailesindeki korunma desenine bakar (yüksek = zararlı).",
    "mutationtaster": "DNA değişiminin hastalık yapma olasılığı (1'e yakın = zararlı).",
    "metalr": "Diğer araçları lojistik regresyonla birleştiren meta skor.",
    "metasvm": "Diğer araçları destek vektör makinesiyle birleştiren meta skor (0 üstü = zararlı).",
    "dann": "CADD'in derin sinir ağı versiyonu.",
    "provean": "Dizi hizalamasıyla değişimin etkisini ölçer (-2,5 altı = zararlı).",
    "vest4": "Hastalık ve zararsız varyantlarla eğitilmiş rastgele orman modeli.",
    "mutpred": "Değişimin hangi mekanizmayı bozduğunu da tahmin eden model.",
    "clinpred": "ClinVar ile eğitilmiş, popülasyon frekansını da kullanan model.",
    "gerp": "DNA bazının memeliler boyunca ne kadar korunduğu (yüksek = önemli yer).",
    "phylop100": "100 omurgalı türünde bu bazın korunmuşluğu (yüksek = önemli yer).",
    "genocanyon": "Genomik bölgenin genel işlevselliği; değişimi değil bölgeyi puanlar.",
    "fitcons": "Bölgenin evrimsel önemi; değişimi değil bölgeyi puanlar.",
    "spliceai": "Değişimin ekson birleşmesini (splicing) bozma olasılığı (0,2 üstü = bozabilir).",
    "dbscsnv_ada": "Ekson sınırına yakın değişimlerde splice bozulması olasılığı.",
    "dbscsnv_rf": "Ekson sınırına yakın değişimlerde splice bozulması olasılığı (rastgele orman).",
}

# Kenar cubugundaki "Bu yazilar nasil okunur?" kutusu
NOTATION = """
Varyantlar iki dille yazılır. **Harf–sayı–harf** proteindeki değişikliği, **c.** ile başlayan yazım DNA'daki
değişikliği anlatır (HGVS kuralı).

**Protein: R175H**
R = normal amino asit (Arginin) · 175 = proteindeki sırası (TP53'te 393 amino asit var) · H = yerine gelen (Histidin).
Uzun yazımı: p.Arg175His.

| Harf | Amino asit | Harf | Amino asit |
|---|---|---|---|
| A | Alanin | M | Metiyonin |
| C | Sistein | N | Asparajin |
| D | Aspartat | P | Prolin |
| E | Glutamat | Q | Glutamin |
| F | Fenilalanin | R | Arginin |
| G | Glisin | S | Serin |
| H | Histidin | T | Treonin |
| I | İzolösin | V | Valin |
| K | Lizin | W | Triptofan |
| L | Lösin | Y | Tirozin |

**Yıldız (\\*): R213\\*** — 213. sırada "dur" sinyali; protein orada biter (nonsense, genelde zararlı).
R213X ve p.Arg213Ter aynı şeydir.

**DNA: c.743G>A** — kodlayan DNA'nın 743. harfi G iken A olmuş. 3 DNA harfi 1 amino asit kodladığı için
743 ÷ 3 ≈ 248 → bu, R248Q'nun DNA yazımıdır.

**c.215del** — 215. harf silinmiş; okuma çerçevesi kayar (frameshift), protein bozulur.

**c.375+1G>A** — 375. harften sonraki ara bölgenin (intron) 1. harfi değişmiş. Ekson sınırındaki ilk 2 harf
"birleştirme" sinyalidir; bozulunca splicing bozulur.

| İşaret | Anlamı | Örnek |
|---|---|---|
| `>` | bir harf başka harfe dönüştü | c.524G>A |
| `del` | silindi | c.215del |
| `dup` | bir kopya fazladan eklendi | c.1179dup |
| `ins` | araya harf eklendi | c.636_637insT |
| `delins` | silinip yerine başka harf kondu | c.743_744delinsAA |
| `fs` | çerçeve kayması | R213fs |
| `=` | amino asit değişmedi (sinonim) | R213= |
| `+5`, `-2` | intronda, eksonun sonundan sonra / başından önce | c.375+5G>A, c.1101-2A>G |
| `c.-29` | proteinin başlangıcından önceki kısım (5' UTR) | c.-29+1G>T |
| `c.*1175` | "dur" sinyalinden sonraki kısım (3' UTR) | c.*1175A>C |
| `_` | aralık: şuradan şuraya | c.637_640del |

**Arama kutusuna yazılabilenler:** R175H · p.Arg175His · c.524G>A · TP53:c.524G>A · chr17-7578406-C-T (hg19) ·
17-7675088-C-T (hg38) · rs28934578 · R213R (sinonim) · M1? (başlangıç kaybı)
"""

# ACMG kanit gruplari: (Turkce baslik, grubun sordugu soru). Anahtarlar Franklin'deki Ingilizce grup adlari.
CATEGORY_TR = {
    "Functional Studies": ("Laboratuvar testleri", "Laboratuvarda protein bozuk mu çıktı, normal mi?"),
    "Segregation Data": ("Aile içi geçiş", "Varyant ailede kanserli kişilerle birlikte mi geçiyor?"),
    "Case Control Studies": ("Hasta sayısı", "Varyant Li-Fraumeni kriterlerine uyan kaç hastada görüldü?"),
    "Population Data": ("Toplumdaki sıklık", "Varyant sağlıklı insanlarda ne kadar sık?"),
    "Effect on Protein": ("Proteine etkisi", "Değişiklik proteini kesiyor mu, bilinen zararlı bir değişime benziyor mu?"),
    "Functional Data": ("Kritik bölge", "Varyant proteinin önemli / sık bozulan bir bölgesinde mi?"),
    "In-silico Predictions": ("Bilgisayar tahminleri", "Tahmin araçları (BayesDel + Align-GVGD) zararlı diyor mu?"),
    "Reputable Source Data": ("Başka kaynakların kararı", "ClinVar gibi kaynaklar ne diyor? (TP53 kurallarında kullanılmaz)"),
    "De Novo Data": ("Yeni ortaya çıkma (de novo)", "Varyant anne-babada yokken çocukta mı ortaya çıktı?"),
    "Allelic Data": ("Diğer kopya", "Aynı kişide başka bir zararlı TP53 varyantı var mı?"),
    "Phenotype": ("Klinik tablo", "Hastanın bulguları bu gene özgü mü? (TP53'te: kanda düşük oranlı varyant)"),
    "Alternate Locus": ("Başka açıklama", "Hastalığı başka bir gen açıklıyor mu?"),
}

# Kriter kodunun sade anlami (ust ozet tablosu icin)
CODE_TR = {
    "PVS1": "Protein kesiliyor / hiç üretilmiyor", "PS1": "Aynı amino asit değişimi zaten patojenik",
    "PS2": "Anne-babada yok, çocukta yeni ortaya çıkmış", "PS3": "Laboratuvarda protein bozuk çıktı",
    "PS4": "Li-Fraumeni hastalarında beklenenden sık", "PM1": "Sık bozulan kritik bölgede (hotspot)",
    "PM2": "Sağlıklı insanlarda çok nadir / hiç yok", "PM4": "Protein uzunluğu değişiyor",
    "PM5": "Aynı pozisyonda başka patojenik değişim var", "PM6": "Varsayılan de novo",
    "PP1": "Ailede kanserle birlikte geçiyor", "PP2": "Bu gende missense genelde patojenik",
    "PP3": "Tahmin araçları zararlı diyor", "PP4": "Kanda düşük oranlı varyant (klonal hematopoez ipucu)",
    "PP5": "Başka kaynak patojenik diyor", "BA1": "Toplumda çok yaygın (tek başına benign)",
    "BS1": "Hastalık için fazla yaygın", "BS2": "60 yaş üstü kansersiz taşıyıcılar var",
    "BS3": "Laboratuvarda protein normal çıktı", "BS4": "Ailede kanserle birlikte geçmiyor",
    "BP1": "Bu gende missense genelde zararsız", "BP2": "Aynı kişide başka zararlı varyant var",
    "BP3": "Tekrar bölgesinde çerçeveyi bozmayan değişim", "BP4": "Tahmin araçları zararsız diyor",
    "BP5": "Hastalığı başka gen açıklıyor", "BP6": "Başka kaynak benign diyor",
    "BP7": "Sinonim ve splice etkisi yok",
}
