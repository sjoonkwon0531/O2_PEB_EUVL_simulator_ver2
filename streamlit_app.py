# -*- coding: utf-8 -*-
"""O2-PEB 공정 탐색기 — 모델 v4 (라디칼 사슬 반응), imec 2026 SPIE 데이터 회귀판.

닫힌형:  (1 - phi_th*x) * exp(-K*x) = 1 - phi_th,   x = D2G / D2G(N2)
         K = k180 * exp[-(Ea/R)(1/T - 1/453.15)] * theta(O2)/theta(21%) * t_O2 * exp(-t_PED,air / tau)
         theta(p) = p / (p + K_O2)
파라미터는 Pollentier et al., Proc. SPIE 13983, 1398312 (2026) 상용 MOR 9점 회귀 (phi_th = 0.5).
"""
import numpy as np
import streamlit as st
import plotly.graph_objects as go

st.set_page_config(page_title="O₂-PEB 공정 탐색기 (모델 v4)", layout="wide")

# ---------------------------------------------------------------- model v4
R = 8.314
LNK0, EA, KO2 = -3.5478346360612534, 113.87321933311503, 58.52905932889255  # fit_com.json, phi_th=0.5
PH, D0, TAU, T0 = 0.5, 60.0, 1440.0, 453.15
AU, EAU = 6.26e2, 71.8                         # Dorney Table 1, unexposed thermolysis
POR = dict(T=180.0, t=60.0, p=21.0)


def x_of(K):
    """Vectorised bisection for (1 - PH x) exp(-K x) = 1 - PH, x in (0, 1]."""
    K = np.atleast_1d(np.asarray(K, dtype=float))
    lo = np.full(K.shape, 1e-12)
    hi = np.ones(K.shape)
    for _ in range(55):
        mid = 0.5 * (lo + hi)
        g = (1 - PH * mid) * np.exp(-K * mid) - (1 - PH)
        pos = g > 0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    x = 0.5 * (lo + hi)
    return np.where(K <= 1e-12, 1.0, x)


def theta(p):
    return p / (p + KO2)


def K_of(T, t, p, f=1.0, ped=0.0):
    return (np.exp(LNK0 - EA * 1e3 / R * (1 / (np.asarray(T) + 273.15) - 1 / T0))
            * theta(np.asarray(p, dtype=float)) * np.asarray(t) * f * np.exp(-ped / TAU))


def d2g(T, t, p, f=1.0, ped=0.0, rh=1.0):
    return D0 * x_of(K_of(T, t, p, f, ped)) * rh


def gain(T, t, p, f=1.0, ped=0.0):
    K = K_of(T, t, p, f, ped)
    return K * x_of(K)


def var_ratio(f, G):
    G = np.atleast_1d(G)
    if f >= 0.999:
        return np.ones_like(G)
    tt = np.linspace(0, 1, 401)[:, None]
    lam_u = G[None, :] * tt
    lam_l = np.where(tt < 1 - f, 0.0, G[None, :] * (tt - (1 - f)) / f)
    vu = np.trapezoid(np.exp(-2 * (G - lam_u)), tt[:, 0], axis=0)
    vl = np.trapezoid(np.exp(-2 * (G - lam_l)), tt[:, 0], axis=0)
    return vl / vu


def ler_change(G, f=1.0):
    """LER_3sigma change (%) at matched CD; half shot noise, half in-process noise."""
    G = np.maximum(np.atleast_1d(G), 1e-9)
    S = G / np.expm1(G)
    V = (1 - np.exp(-2 * G)) / (2 * G)
    r2 = 0.5 * V / S ** 2 * var_ratio(f, G) + 0.5 * (np.sinh(G / 2) / (G / 2)) ** 2
    return (np.sqrt(r2) - 1) * 100


def sens(T, t, p, f=1.0, ped=0.0, sT=0.2, sp=0.3, st_=1.0):
    """Dose-equivalent noise (%) from hotplate T, O2 and time control errors."""
    l0 = np.log(d2g(T, t, p, f, ped))
    gT = (np.log(d2g(np.asarray(T) + 0.1, t, p, f, ped)) - l0) / 0.1
    gp = np.where(np.asarray(p) > 0, (np.log(d2g(T, t, np.asarray(p) + 0.1, f, ped)) - l0) / 0.1, 0.0)
    gt = (np.log(d2g(T, np.asarray(t) + 0.5, p, f, ped)) - l0) / 0.5
    return np.sqrt((gT * sT) ** 2 + (gp * sp) ** 2 + (gt * st_) ** 2) * 100


def dark(T, t):
    return (1 - np.exp(-AU * np.exp(-EAU * 1e3 / (R * (np.asarray(T) + 273.15))) * np.asarray(t))) * 100


POR_D = float(d2g(**POR)[0])
POR_S = float(sens(POR["T"], POR["t"], POR["p"])[0])

# ---------------------------------------------------------------- UI
st.title("O₂-PEB 공정 탐색기 — 모델 v4")
st.caption("라디칼 사슬 반응 모델을 imec 2026 SPIE 데이터(상용 MOR)에 회귀한 결과로 PEB 온도·시간·O₂ 분압·"
           "O₂ 공급 시점·PED 분위기를 바꿔 봅니다. 값은 상대 트렌드이며, 라인 소재로는 하이닉스 스플릿 데이터로 "
           "재보정해야 합니다.")

if "T" not in st.session_state:
    st.session_state.update(T=175, t=110, p=21, f=1.0, ped=0)


def set_por():
    st.session_state.update(T=180, t=60, p=21, f=1.0, ped=0, rh="약 5%RH (imec 회귀 조건)")


with st.sidebar:
    st.header("공정 조건")
    T = st.slider("PEB 온도 (°C)", 160, 230, key="T")
    t = st.slider("PEB 시간 (s)", 30, 180, step=5, key="t")
    p = st.slider("PEB O₂ 농도 (%)", 0, 60, key="p",
                  help="공기 = 21%. 23.5% 초과는 산소 농화 분위기(OSHA 29 CFR 1910.146)로 별도 안전 설비가 필요합니다.")
    f = st.slider("O₂ 공급 구간 (베이크 끝 쪽 비율)", 0.2, 1.0, step=0.05, key="f",
                  help="1.0 = 베이크 전체 균일 공급. 작을수록 베이크 후반에만 O₂를 넣습니다.")
    ped = st.slider("공기 중 PED 시간 (s)", 0, 600, step=10, key="ped",
                    help="0 = 진공·N₂ 중 대기. 공기 중 대기는 사슬 운반체를 잃게 합니다 (τ ≈ 24분).")
    rh_lab = st.selectbox("PEB 습도", ["약 5%RH (imec 회귀 조건)", "약 45%RH (상용 MOR 실측 약 −10%)"], key="rh")
    rh = 1.0 if rh_lab.startswith("약 5") else 0.9
    st.header("제약 조건")
    ler_lim = st.slider("LER 허용 증가 (%)", 1.0, 10.0, 3.0, 0.5)
    rob = st.slider("dose 산포 허용 (현행 대비 배수)", 1.0, 2.0, 1.0, 0.05,
                    help="산포 = 온도 ±0.2 °C, O₂ ±0.3%p, 시간 ±1 s에서 생기는 dose 변동")
    tmax = st.slider("PEB 시간 상한 (s)", 60, 180, 120, 5)
    st.button("현행 조건으로 (180 °C, 60 s, 공기)", on_click=set_por)

# ---------------------------------------------------------------- current point
D = float(d2g(T, t, p, f, ped, rh)[0])
G = float(gain(T, t, p, f, ped)[0])
LE = float(ler_change(G, f)[0])
SE = float(sens(T, t, p, f, ped)[0])
DK = float(dark(T, t))

c = st.columns(6)
c[0].metric("D2G (겔 dose)", f"{D:.1f} mJ/cm²", f"{(D / POR_D - 1) * 100:+.0f}% vs 현행", delta_color="inverse")
c[1].metric("루프 게인 G", f"{G:.2f}", "O₂ 사슬 기여", delta_color="off")
c[2].metric("LER 변화 (같은 CD)", f"{LE:+.1f}%", "충족" if LE <= ler_lim else "초과",
            delta_color="normal" if LE <= ler_lim else "inverse")
c[3].metric("dose 산포", f"{SE:.2f}%", f"현행 {POR_S:.2f}%", delta_color="off")
c[4].metric("비노광부 리간드 손실", f"{DK:.3f}%", "Dorney 비노광 속도식", delta_color="off")
c[5].metric("O₂ 분위기", f"{p}%", "산소 농화" if p > 23.5 else "대기 수준", delta_color="off")


# ---------------------------------------------------------------- optimiser
@st.cache_data(show_spinner=False)
def optimise(p_list, f, ped, rh, ler_lim, rob, tmax):
    best = None
    Tg = np.arange(160, 231, 1.0)
    tg = np.arange(30, tmax + 1, 5.0)
    TT, tt = np.meshgrid(Tg, tg)
    for pp in p_list:
        Gm = gain(TT.ravel(), tt.ravel(), pp, f, ped)
        ok = (ler_change(Gm, f) <= ler_lim) & (sens(TT.ravel(), tt.ravel(), pp, f, ped) <= POR_S * rob + 1e-9)
        if not ok.any():
            continue
        Dm = d2g(TT.ravel(), tt.ravel(), pp, f, ped, rh)
        k = np.argmin(np.where(ok, Dm, np.inf))
        if best is None or Dm[k] < best[0]:
            best = (float(Dm[k]), float(TT.ravel()[k]), float(tt.ravel()[k]), float(pp))
    return best


o1, o2 = st.columns(2)
best_p = optimise((float(p),), f, ped, rh, ler_lim, rob, tmax)
best_all = optimise((0.0, 10.0, 15.0, 21.0, 23.5, 30.0, 35.0, 40.0, 45.0, 50.0), f, ped, rh, ler_lim, rob, tmax)
with o1:
    if best_p:
        st.success(f"**현재 O₂ {p}%에서 최적점:** {best_p[1]:.0f} °C, {best_p[2]:.0f} s → "
                   f"D2G {best_p[0]:.1f} mJ/cm² (현행 대비 {(best_p[0] / POR_D - 1) * 100:+.0f}%)")
    else:
        st.warning("현재 O₂에서 제약을 모두 만족하는 조건이 없습니다. 허용폭을 넓혀 보세요.")
with o2:
    if best_all:
        st.info(f"**O₂까지 함께 최적화:** {best_all[1]:.0f} °C, {best_all[2]:.0f} s, O₂ {best_all[3]:.0f}% → "
                f"D2G {best_all[0]:.1f} mJ/cm² (현행 대비 {(best_all[0] / POR_D - 1) * 100:+.0f}%)"
                + (" · 산소 농화 설비 필요" if best_all[3] > 23.5 else ""))

# ---------------------------------------------------------------- maps
@st.cache_data(show_spinner=False)
def map_Tt(p, f, ped, rh):
    Tg = np.linspace(160, 230, 57)
    tg = np.linspace(30, 180, 51)
    TT, tt = np.meshgrid(Tg, tg)
    Dm = d2g(TT.ravel(), tt.ravel(), p, f, ped, rh).reshape(TT.shape)
    Lm = ler_change(gain(TT.ravel(), tt.ravel(), p, f, ped), f).reshape(TT.shape)
    Sm = sens(TT.ravel(), tt.ravel(), p, f, ped).reshape(TT.shape)
    return Tg, tg, Dm, Lm, Sm


@st.cache_data(show_spinner=False)
def map_pT(t, f, ped, rh):
    pg = np.linspace(0, 60, 49)
    Tg = np.linspace(160, 230, 57)
    PP, TT = np.meshgrid(pg, Tg)
    Dm = d2g(TT.ravel(), t, PP.ravel(), f, ped, rh).reshape(PP.shape)
    Lm = ler_change(gain(TT.ravel(), t, PP.ravel(), f, ped), f).reshape(PP.shape)
    Sm = sens(TT.ravel(), t, PP.ravel(), f, ped).reshape(PP.shape)
    return pg, Tg, Dm, Lm, Sm


def heat(xg, yg, Dm, Lm, Sm, xlab, ylab, cur, opt, por, hline=None):
    fig = go.Figure()
    fig.add_heatmap(x=xg, y=yg, z=Dm, colorscale="Viridis", zmin=0, zmax=float(np.nanmax(Dm)),
                    colorbar=dict(title="D2G<br>mJ/cm²"),
                    hovertemplate=f"{xlab}: %{{x:.0f}}<br>{ylab}: %{{y:.0f}}<br>D2G: %{{z:.1f}}<extra></extra>")
    fig.add_contour(x=xg, y=yg, z=Dm, showscale=False, hoverinfo="skip",
                    contours=dict(coloring="none", showlabels=True, start=5, end=60, size=5,
                                  labelfont=dict(color="white", size=10)),
                    line=dict(color="rgba(255,255,255,0.6)", width=1))
    fig.add_contour(x=xg, y=yg, z=Lm, showscale=False, hoverinfo="skip", name="LER 한계",
                    contours=dict(coloring="none", start=ler_lim, end=ler_lim, size=1),
                    line=dict(color="#ff5c5c", width=3))
    fig.add_contour(x=xg, y=yg, z=Sm, showscale=False, hoverinfo="skip", name="산포 한계",
                    contours=dict(coloring="none", start=POR_S * rob, end=POR_S * rob, size=1),
                    line=dict(color="#ffd84d", width=3, dash="dash"))
    if hline:
        fig.add_hline(y=hline, line=dict(color="white", dash="dot", width=1.5))
    fig.add_scatter(x=[por[0]], y=[por[1]], mode="markers", name="현행",
                    marker=dict(symbol="diamond", size=13, color="white", line=dict(color="black", width=1.5)))
    if opt:
        fig.add_scatter(x=[opt[0]], y=[opt[1]], mode="markers", name="최적점",
                        marker=dict(symbol="star", size=18, color="#ff2fd0", line=dict(color="black", width=1.5)))
    fig.add_scatter(x=[cur[0]], y=[cur[1]], mode="markers", name="현재",
                    marker=dict(symbol="circle", size=12, color="white", line=dict(color="black", width=2)))
    fig.update_layout(height=470, margin=dict(l=10, r=10, t=10, b=10), xaxis_title=xlab, yaxis_title=ylab,
                      legend=dict(orientation="h", y=-0.18))
    return fig


m1, m2 = st.columns(2)
with m1:
    st.subheader("D2G 지도: 온도 × 시간")
    st.caption("현재 O₂·공급 구간·PED 기준. 빨강 = LER 한계, 노랑 점선 = 산포 한계, 흰 점선 = 시간 상한")
    Tg, tg, Dm, Lm, Sm = map_Tt(float(p), f, float(ped), rh)
    st.plotly_chart(heat(Tg, tg, Dm, Lm, Sm, "PEB 온도 (°C)", "PEB 시간 (s)", (T, t),
                         (best_p[1], best_p[2]) if best_p else None, (180, 60), hline=tmax),
                    width="stretch")
with m2:
    st.subheader("D2G 지도: O₂ × 온도")
    st.caption("현재 PEB 시간 기준. O₂ 농도와 온도의 등가 관계가 보입니다.")
    pg, Tg2, Dm2, Lm2, Sm2 = map_pT(float(t), f, float(ped), rh)
    st.plotly_chart(heat(pg, Tg2, Dm2, Lm2, Sm2, "PEB O₂ 농도 (%)", "PEB 온도 (°C)", (p, T),
                         (best_all[3], best_all[1]) if best_all and best_all[2] == t else None, (21, 180)),
                    width="stretch")

# ---------------------------------------------------------------- model notes
with st.expander("모델 v4 — 방정식과 파라미터 (공정 단위)", expanded=False):
    st.latex(r"(1-\varphi_{th}\,x)\,e^{-Kx} = 1-\varphi_{th},\qquad x=\mathrm{D2G}/\mathrm{D2G}(\mathrm{N_2})")
    st.latex(r"K = k_{180}\,e^{-\frac{E_a}{R}\left(\frac{1}{T}-\frac{1}{453\,\mathrm{K}}\right)}\,"
             r"\frac{\theta(\mathrm{O_2})}{\theta(21\%)}\,t_{\mathrm{O_2}}\,e^{-t_{\mathrm{PED,air}}/\tau},\qquad "
             r"\theta=\frac{\mathrm{O_2}}{\mathrm{O_2}+K_{\mathrm{O_2}}}")
    st.markdown(f"""
| 기호 | 뜻 | 값 |
|---|---|---|
| D2G(N₂) | O₂ 없이 겔이 되는 노광량 | {D0:.0f} mJ/cm² |
| φ_th | 겔 문턱 전환율 | {PH} (가정, 하이닉스 contrast curve로 고정 예정) |
| k₁₈₀ | 180 °C·공기의 사슬 반응 속도 | {np.exp(LNK0) * 21 / (21 + KO2):.4f} s⁻¹ |
| Eₐ | 사슬 전파 활성화에너지 | {EA:.0f} ± 10 kJ/mol (+10 °C마다 ×1.9) |
| K_O₂ | O₂ 효과 반포화 농도 | {KO2:.0f} ± 42 % |
| τ | 공기 중 PED 운반체 손실 시간 | {TAU:.0f} s |

**기전:** EUV가 Sn–C를 끊어 생긴 Sn 라디칼을 O₂가 장벽 없이 포획(Bu₃Sn• + O₂ = 7.5×10⁹ M⁻¹s⁻¹,
Maillard·Ingold·Scaiano, JACS 1983)하고, 그 peroxyl이 이웃 리간드를 공격하는 사슬 반응(Eₐ ≈ 59 kJ/mol, FTIR).
측정되는 빈도인자는 A_p × θ(O₂) × [Sn 라디칼]이라 O₂는 Eₐ를 바꾸지 않고 A만 키우며, 비노광부에는 영향이 없습니다.

**LER:** 같은 CD 닫힌형(shot·공정 노이즈 반반). O₂ 공급 구간은 공정 노이즈 감쇠 시점 효과(이론 노트 식 8)로 반영.

**데이터:** Pollentier et al., Proc. SPIE 13983, 1398312 (2026); Dorney et al., Proc. SPIE 13983, 139831G (2026).
""")
