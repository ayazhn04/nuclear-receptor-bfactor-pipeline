# SASA + IQR: результаты для команды (29 октября 2026)

**Новый exploratory-анализ исправленного v2; старые результаты не заменяет.**

SASA и B-factor сопоставлены для 416,800 PDB × residue наблюдений. Основной порог: ≥5 PDB с полными тяжёлыми атомами; доступны 7,065 позиций в 31 рецепторах.

Координаты: 1771/1776 PDB совпали с frozen SHA256. Без verified coordinates: 29OB, 6PIT, 9W11, 9W12, 9W13. Calculation failures: нет. Эти записи сохранены в master, но не участвуют в новом SASA-анализе.

## 1. Связь экспозиции с B-factor и IQR

Корреляция Spearman рассчитана отдельно в каждом рецепторе (≥20 позиций). Доверительный интервал — bootstrap рецепторов, не остатков. Wilcoxon проверяет распределение receptor-level rho относительно нуля; BH для двух endpoints.

- IQR_Z: n=30, медианная rho=0.414, 95% CI [0.368; 0.459], q=1.863e-09.
- MEDIAN_Z: n=30, медианная rho=0.560, 95% CI [0.536; 0.576], q=1.863e-09.

## 2. Лигандный карман: IQR после matching

Сравниваются LBP и позиции вне LBP/AF2-H12 одного рецептора: точный professor-контекст, support 5–9/10–19/≥20, затем RSA ≤0.05 / (0.05,0.25] / >0.25. ≥1 site и ≥3 controls на страту; медиана разностей по стратам → один эффект на рецептор. Отрицательный эффект = меньший IQR у LBP. P — receptor-level Wilcoxon; BH для четырёх designs.

- **CONTEXT_SUPPORT**: n=29; site positions=750; ΔIQR=-0.0329, 95% CI [-0.0654; 0.0018]; p=0.0592; q=0.2368.
- **CONTEXT_SUPPORT_RSA**: n=29; site positions=604; ΔIQR=0.0020, 95% CI [-0.0523; 0.0333]; p=0.7983; q=1.
- **CONTEXT_SUPPORT_RSA_NO_REVIEW**: n=29; site positions=591; ΔIQR=0.0020, 95% CI [-0.0523; 0.0380]; p=0.848; q=1.
- **CONTEXT_SUPPORT_RSA_AMINOACID**: n=6; site positions=8; ΔIQR=-0.0146, 95% CI [-0.0278; 0.0958]; p=1; q=1.

## Вывод для презентации

После учёта экспозиции и структурного контекста убедительного LBP-специфичного различия IQR не установлено по совместному критерию CI вне нуля и BH q<0.05. Это не доказывает отсутствия эффекта.

**IQR(z) ≠ median(z).** IQR показывает неоднородность относительного кристаллографического disorder между PDB-записями, а не абсолютную гибкость или динамику в растворе. Корреляции с median(z) и с IQR(z) отвечают на разные вопросы.

## Что нельзя утверждать

- Нельзя объявлять ligand-induced эффект: это cross-sectional association, не APO/HOLO intervention.
- AF2/H12 описан по IQR и RSA, но exact-context site-specific тест не выполнялся: весь H12 определяет сайт, нет независимого same-H12 background.
- Это isolated-chain SASA, не SASA комплекса/биологической сборки/кристаллической решётки; ligand/partner occlusion не измерена.
- Неполные side chains и неканонические residue types не получают RSA. Missingness и coarse bins могут влиять на результат.
- Разрешение, refinement, packing, ligand composition и пространственная зависимость позиций полностью не контролируются; confidence intervals не снимают эти ограничения.
- Matching меняет набор доступных рецепторов/позиций; сравнение designs не является формальной оценкой mediation/attenuation.

## Графики

1. `01_sasa_correlations.png` — receptor-level корреляции.
2. `02_sasa_adjusted_lbp_iqr.png` — эффекты matching и CIs.
3. `03_sasa_iqr_examples.png` — шесть рецепторов с максимальным числом eligible positions, без pooled-тестов.
4. `04_functional_site_iqr_rsa_heatmap.png` — IQR/RSA функциональных сайтов по рецепторам, описательно.

Методы, команды и дополнительные ограничения: `README.md`. Provenance: `sasa_manifest.json`, `analysis_manifest.json` и audit CSV.
