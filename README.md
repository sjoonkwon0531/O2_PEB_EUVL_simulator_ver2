# O₂-PEB 공정 탐색기 (모델 v4)

High NA EUV 금속 산화물 레지스트(MOR)의 O₂-assisted PEB를 **라디칼 사슬 반응 모델(v4)**로 계산하는 Streamlit 앱입니다.
모델 파라미터는 imec 2026 SPIE 데이터(Pollentier et al., Proc. SPIE 13983, 1398312; Dorney et al., Proc. SPIE 13983, 139831G)에
회귀했습니다. SK하이닉스 산학과제 "High NA EUV MOR/UL 멀티스케일 전산 설계" (성균관대 SPMDL).

| 페이지 | 내용 |
|---|---|
| 메인 (`streamlit_app.py`) | 공정 탐색기: PEB 온도·시간·O₂ 농도·O₂ 공급 구간·공기 중 PED·습도와 제약(LER, dose 산포, 시간 상한)을 바꾸면 D2G·루프 게인·LER 변화·dose 산포·비노광 손실·O₂ 안전 등급, 온도×시간 / O₂×온도 지도, 제약 하 최적점이 즉시 계산됨 |
| `pages/1_model_notes.py` | imec 핵심 데이터, 9월 모델 대비 수정 사항, 반응 기전과 빈도인자 해석, 회귀 결과, 최적화 결과, 한계와 참고문헌 |

## 모델 요약

```
(1 − φ_th·x)·exp(−K·x) = 1 − φ_th,     x = D2G / D2G(N₂)
K = k₁₈₀ · exp[−(Eₐ/R)(1/T − 1/453 K)] · θ(O₂)/θ(21%) · t_O₂ · exp(−t_PED,air / τ)
θ(O₂) = O₂ / (O₂ + K_O₂)
```

상용 MOR 회귀값: D2G(N₂) = 60 mJ/cm², Eₐ = 114 ± 10 kJ/mol, K_O₂ = 59 ± 42 %, k₁₈₀ = 0.0076 s⁻¹, τ = 1,440 s (φ_th = 0.5 가정).
LER 변화는 같은 CD 닫힌형(shot·공정 노이즈 반반), O₂ 공급 구간은 공정 노이즈 감쇠 시점 효과로 반영합니다.

## 로컬 실행

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

## Streamlit Community Cloud 배포

1. 이 저장소를 GitHub에 올립니다 (private 가능).
2. https://share.streamlit.io → Create app → 저장소/브랜치 선택, Main file path에 `streamlit_app.py`.

## 저장소 구조

```
o2peb-v4-explorer/
├── streamlit_app.py          # 공정 탐색기 (모델 v4)
├── pages/
│   └── 1_model_notes.py       # 근거와 해설 (모델·데이터·최적화)
├── assets/                   # 그림 6장 (imec 데이터, 수정 사항, 기전, 회귀, 최적화)
├── requirements.txt
└── README.md
```

## 주의

연구용 트렌드 도구입니다. 모든 값은 imec 상용 MOR 기준 상대값이며, 라인 소재 적용 전 하이닉스 스플릿 데이터(contrast curve,
O₂·온도·시간 스플릿)로 재보정해야 합니다.
