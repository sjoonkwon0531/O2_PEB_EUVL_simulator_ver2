# -*- coding: utf-8 -*-
"""모델 v4의 근거: imec 데이터, 9월 모델 대비 수정 사항, 반응 기전, 회귀, 최적화 결과."""
import os
import streamlit as st

st.set_page_config(page_title="모델 v4 — 근거와 해설", layout="wide")
A = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets")


def img(name, caption):
    st.image(os.path.join(A, name), caption=caption, width="stretch")


st.title("모델 v4 — 근거와 해설")
st.caption("imec 2026 SPIE 두 편(Pollentier 13983-36: 대비곡선 D2G, Dorney 13983-50: in-situ FTIR)으로 "
           "9월 모델(v2)을 고치고, 라디칼 사슬 반응 기전 모델(v4)로 올린 과정입니다.")

st.header("1. imec 핵심 데이터")
img("FI6_imec_data.png",
    "a: 상용 MOR, PED/PEB 분위기 분리 (Pollentier 표 1, 180 °C, 각 120 s) — PEB 분위기가 결과를 결정. "
    "b: Dorney 표 1의 A·Eₐ로 재계산한 리간드 소모 속도 — 노광부만 O₂로 평행 이동, 비노광부는 O₂ 무관")
st.markdown("""
| 관찰 | 수치 |
|---|---|
| O₂ 21→50% | D2G −15% (상용 MOR, 180 °C), −25% (BMOR), −30% (OSMO) |
| O₂ 유무 비교 | O₂를 빼면 상용 MOR D2G 약 2배, 모델 MOR은 불용화 없음 |
| O₂ 없을 때 온도·시간 | 160→200 °C, 60→240 s에도 D2G 거의 불변 (O₂ 있으면 −55~75%) |
| Arrhenius | 노광부 Eₐ 58.6 kJ/mol (21%·50% 동일), A만 증가. 비노광부 O₂ 무관 |
| FTIR | C=O(1,700 cm⁻¹)는 O₂·온도에 비례, CO₂·H₂O 무관 |
""")

st.header("2. 9월 모델(v2)의 수정 사항")
img("FV1_old_vs_new.png",
    "a: O₂ 의존성 — v2는 N₂에서 dose +15% 예측, 실측 +69%. b: PEB 시간 — v2는 O₂ 없이도 시간 4배면 dose 1/4 예측, "
    "실측 0.90. c: 활성화에너지 — v2 예측 30–38 kJ/mol은 낮았음, 라디칼 소거(Scenario A)는 기각")
st.markdown("""
| 9월 가정 / 예측 | 판정 | 새 모델 (v4) |
|---|---|---|
| O₂ = 라디칼 소거 → dose 증가 (Scenario A) | **폐기** | O₂ 포획은 사슬을 시작하는 "생산적 포획" |
| Eₐ 30–38 kJ/mol | **폐기** | 사슬 전파 단계 Eₐ 59 (FTIR) / 114 (D2G 회귀) |
| 베이크 전 흡수(pre-soak) 이득 | **폐기** | 상온에서는 전파가 거의 멈춤 → PEB 중 공급이 필요 |
| N₂ 대비 공기 이득 +15–23% | **폐기** | 실측 약 +70%: O₂ 경로가 주 경로 |
| 기저 반응이 PEB 시간에 비례 | 수정 | 기저 = 즉시 완료 (EUV 직접 절단) |
| ηc ∝ O₂ 선형 | 수정 | θ = O₂/(O₂ + K_O₂), K_O₂는 소재 의존 |
| 피드백 a(x)·ηc(1−φ), 용량 한계, LER 불변, 꼬리 편집 | 유지 | 그대로 (v4에서 재유도 예정) |
""")

st.header("3. 반응 기전과 빈도인자")
img("FV2_mechanism.png",
    "a: 반응 좌표 개략도 — O₂ 포획은 장벽 없음, 사슬 전파에 59 kJ/mol 장벽. "
    "b: O₂ 10·21·50%의 Arrhenius 직선이 같은 기울기로 평행 이동 = 'A만 증가'")
st.latex(r"k_{obs} = k_p\,[\mathrm{Sn\!-\!OO^\bullet}] = A_p\,e^{-E_a/RT}\cdot\theta(\mathrm{O_2})\cdot s,\qquad "
         r"\theta=\frac{\mathrm{O_2}}{\mathrm{O_2}+K_{\mathrm{O_2}}}")
st.markdown("""
- **왜 Eₐ는 같고 A만 커지나:** O₂는 장벽 없는 포획 단계에만 관여합니다(Bu₃Sn• + O₂ = 7.5×10⁹ M⁻¹s⁻¹,
  Maillard·Ingold·Scaiano, *J. Am. Chem. Soc.* 105 (1983), doi:10.1021/ja00353a039). 장벽은 그 뒤 사슬 전파 단계에 있으므로,
  O₂는 언덕 높이가 아니라 언덕을 넘으려는 운반체의 수를 바꿉니다.
- **왜 노광부만:** 운반체(Sn 라디칼)는 EUV가 Sn–C를 끊어야 생깁니다. 비노광부에는 포획할 상대가 없습니다.
- **왜 O₂에 비례하지 않나:** 포획이 라디칼 손실과 경쟁해 포화합니다. 상용 MOR K_O₂ ≈ 59% (±42), OSMO는 선형 — 소재 의존.
- **Dorney의 A 비 1.46:** Eₐ 오차(±13 kJ/mol)와 강하게 상관되어(A의 불확도 약 27배) O₂ 차수를 정하는 데 부적합합니다.
""")

st.header("4. imec 데이터 3종 동시 재현")
img("FV3_joint_fit.png",
    "a: 상용 MOR 4개 온도·시간 조건 (rms 7.8%). b: OSMO (θ 선형 영역). c: Dorney 고정 온도 O₂ 스캔과 같은 함수형. "
    "d: PED/PEB 표 — vac/air는 예측값 (25 대 실측 22)")
st.markdown("""
| 기호 | 뜻 | 값 (상용 MOR, φ_th = 0.5) |
|---|---|---|
| D2G(N₂) | O₂ 없이 겔이 되는 노광량 | 60 mJ/cm² |
| k₁₈₀ | 180 °C·공기의 사슬 반응 속도 | 0.0076 s⁻¹ (60 s 베이크 = 루프 게인 약 0.3) |
| Eₐ | 사슬 전파 활성화에너지 | 114 ± 10 kJ/mol (+10 °C마다 ×1.9) |
| K_O₂ | O₂ 효과 반포화 농도 | 59 ± 42 % |
| τ | 공기 중 PED 운반체 손실 시간 | 1,440 s |

등가 관계: O₂ 21→50% ≈ PEB 온도 +8 °C ≈ PEB 시간 ×1.7 (같은 dose).
""")

st.header("5. 공정 최적화 결과")
img("FV4_coopt_maps.png",
    "O₂ 21%(위)·50%(아래)의 온도×시간 지도. a·d: D2G, b·e: 같은 CD의 LER 변화, c·f: dose 산포. "
    "★ = LER +3%·현행 강건성·120 s 이내의 최소 dose 조건")
img("FV5_profile_ped.png",
    "a: 허용 산포에 따른 최소 D2G (강건성–dose 트레이드오프). b: 베이크 후반 O₂ 공급의 공정 노이즈 감소. "
    "c: 공기 중 PED 대기 시간에 따른 dose 손해")
st.markdown("""
| 레시피 | PEB | O₂ / 공급 | PED | D2G | dose 산포 |
|---|---|---|---|---|---|
| 현행 | 180 °C, 60 s | 21% / 균일 | 공기 | 39.2 | 0.91% |
| A. 저온·장시간 | 175 °C, 112 s | 21% / 균일 | N₂·진공 | 34.6 (−12%) | 0.90% |
| B. 산소 농화 | 170 °C, 118 s | 50% / 균일 | N₂·진공 | 30.5 (−22%) | 0.90% |
| C. 공기, 산포 ×1.5 | 189 °C, 120 s | 21% / 균일 | N₂·진공 | 19.2 (−51%) | 1.36% |
| B′. 후반 공급 | 175 °C, 110 s | 50% / 끝 40% | N₂·진공 | 40.5 (+3%) | 공정 노이즈 σ² −15% |

강건성: 온도 ±0.2 °C, O₂ ±0.3%p, 시간 ±1 s에서 생기는 dose 변동(가정값). 장시간 쪽이 유리한 이유는 시간 오차의 상대 비중이
줄어들기 때문이며, 고온으로 dose를 줄이면 세 항이 모두 커집니다.
""")

st.header("6. 한계와 참고문헌")
st.markdown("""
- 모든 값은 imec 상용 MOR 기준의 상대 트렌드입니다. 라인 소재로는 하이닉스 스플릿 데이터로 재보정이 필요합니다.
- φ_th는 가정값(0.5)이며 contrast curve로 고정해야 합니다. 160 °C에서는 이득을 약간 과대 예측합니다.
- 습도는 경험식(상용 MOR 45%RH 약 −10%), C=O 부반응의 현상 영향은 미반영, imec 데이터 범위(160–220 °C, O₂ ≤ 50%) 밖은 외삽입니다.
- 산소 농화 분위기 기준: OSHA 29 CFR 1910.146 (23.5% 초과).

**참고문헌**
- I. Pollentier et al., Proc. SPIE 13983, 1398312 (2026), doi:10.1117/12.3091787
- K. M. Dorney et al., Proc. SPIE 13983, 139831G (2026), doi:10.1117/12.3091713
- B. Maillard, K. U. Ingold, J. C. Scaiano, J. Am. Chem. Soc. 105(15) (1983), doi:10.1021/ja00353a039
- "Absolute rate constants for hydrocarbon autoxidation. XXI," Can. J. Chem. (1972), doi:10.1139/v72-365
""")
