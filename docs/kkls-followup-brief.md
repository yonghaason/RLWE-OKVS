# KKLS 후속 연구 — 세션 핸드오프 브리프

> **[이관 주석]** 이 브리프와 후속 note 작업은 Real-PJC 레포에서 이 레포(RLWE-OKVS,
> KKLS 구현 레포)로 옮겨졌다. 아래 본문이 참조하는 자료 경로는 여전히
> **yonghaason/Real-PJC 레포 기준**이다: `docs/encrypted-alignment-design.tex`,
> `docs/realpjc.tex`의 git 히스토리. 또한 본문에 적힌 삭제 커밋 해시 `2d828a1`은
> stale 해시로, 실제 삭제 커밋은 Real-PJC의 `043677f`이다
> (`git show 043677f^:docs/realpjc.tex`로 복원). 구현 디렉토리는 이 레포의
> `experiments/`, `rlwe-okvs/`, `exec/` (본문의 `real_okvs/`는 옛 이름).
> 회수 결과물은 `docs/kkls-followup-note.tex`에 이미 반영되어 있다.

> 이 문서는 Real-PJC 논문 세션에서 분기된 핸드오프입니다. **Real-PJC 세션은 Private Data
> Alignment 논문에 집중**하고, 이 브리프를 읽는 새 세션은 **KKLS(base RPMT 논문,
> `EPRINT:KwaKwoLeeSon26`)의 후속 short paper/note**를 담당합니다. 두 주제를 다룹니다:
> (A) BGV의 2×(N/2) slot 구조를 온전히 활용하는 dense two-row decode,
> (B) ssPMT의 layout 공개 최적화에 대한 정확한 leakage 분석과 보수.
>
> 아래 내용은 원 세션에서 도출·검증된 결론의 요약이며, "검증 필요" 표시가 붙은 항목은
> 아직 독립 검증이 안 된 세션 유도입니다.

## 공통 배경 (KKLS 요약)

- Batch homomorphic OKVS decoding: OKVS decode는 linear이므로 배치 decode는 plain
  matrix × ciphertext vector. RSB(spaced-band) 컬럼 순열 π_H로 각 band의 w개 계수를
  연속 블록에 하나씩·공통 intra-block slot에 배치 → 곱이 slot-aligned PMult/Add만으로
  수행됨 (wrap 경계 제외).
- Sequencing: band start가 같은 슬롯을 요구하는 query 충돌을 L개 layer로 분할해 해소.
  query j ↦ (μ_j, τ_j) (layer, slot). 배정은 H₁(F(y_j)), slot 수, span 한계 ζ에만 의존.
  잘 채워진 영역(H ≪ n_y)에서 L ≈ n_y/H.
- ssPMT (KKLS Appendix D 스타일): S가 상수 indicator IS로 OKVS 인코딩, R이 자기 순서로
  homomorphic decode, 마스킹 후 반환, S 복호 → ss-equality → XOR share.
  **핵심 성질: 암호화된 decode라 상수 target이 허용되어, share가 R의 decode layout에서
  "태어나면서부터 정렬"됨** (평문 OPPRF에서는 상수 target이 R에게 반복값으로 새므로 불가).

## Topic A: 2×(N/2) dense two-row decode

### 자료 위치
- **1차 자료: `docs/encrypted-alignment-design.tex`** — 다음이 온전히 남아 있음:
  - "Pre-rotation vs. homomorphic rotation, and the BGV two-key dense layout" 문단
  - "The BGV dense two-row decode" 문단 + Algorithm `alg:homdecode-bgv`
  - "Output-swap refinement" 문단 (rotation 수 L/2+(w−1)/2)
  - Open problem "2-row decode — resolved" (integrated sequencing 잔여 과제 포함)
- **제거된 논문 텍스트**: `docs/realpjc.tex`의 git 히스토리. 커밋 `2d828a1`이 §3.1
  "Variant with Receiver-side Rotation" + Algorithm + "Optimality of the rotation count" +
  Appendix "Redundancy Sweep"(ε-sweep 실측 테이블)을 삭제함.
  복원: `git show 2d828a1^:docs/realpjc.tex`.
- 구현: 레포의 `experiments/`, `real_okvs/`, `exec/` (BGV indicator decode 두 layout 모두
  구현·실측 완료 상태라고 기록됨).

### 기술 요약
- BGV의 N개 slot은 Z_{N/2}×Z_2 hypercube (row rotation + column swap). 평평한 N-cyclic
  shift는 단일 automorphism이 아님 → 단순 1-key 접근은 N/2 slot만 쓰거나 낭비.
- Dense two-row: 입력 ct당 두 개의 H-block 패킹 C_t = Enc(p_{2t} ∥ p_{2t+1}), 출력 ct당
  두 sequenced layer(짝수 layer→row 0, 홀수→row 1). Band가 양쪽 parity 블록에 걸치는
  문제를 column swap σ(C_t)로 라우팅.
  - Input-swap 스케줄: σ(C_t)를 입력당 1회 계산·재사용 → m/N + (w−1) rotations.
  - **Output-swap 스케줄 (기본값)**: 항등식 D⊙σ(x) = σ(σ(D)⊙x)로 swap을 출력측으로
    접어 L/2 + (w−1)/2 rotations, key switch가 plain-multiply 합 뒤에 떨어져 noise
    budget ~1비트 이득.
- 실측 (n=2^20, N=2^13): pre-rotated 1.8s/53MB, dense(output-swap) 1.8s/49MB, base RPMT
  reference 2.0s/58MB. ε-sweep: pre-rotated 통신은 ε≈0.4에서 U자 바닥(43MB), dense는
  단조 감소(ε=0.2에서 28.5MB, −38%); 시간 교차점 ε≈0.6–0.7.
- 잔여 과제: **integrated sequencing** — 현재 구현은 N/2에서 sequence 후 pairing
  (sequence-then-pair). 두 row를 하나의 공유 블록 윈도우(≤ζ)에서 같이 채우면 PMult 수
  0.97–1.04× (vs 현행 1.18–1.41×)로 줄일 수 있음 (측정 기반 추정).

### 후속 논문이 할 일
1. Dense two-row 구성의 정식 서술 + 정확성 증명 (row-wise, merge 없음)
2. Rotation-schedule 공간 정리와 optimality 논증 (input-swap vs output-swap vs
   replicated-p vs per-output merge; L/2 < m/N+(w−1)/2 조건)
3. Integrated sequencing 설계·구현·재실측
4. ε-sweep 테이블 복원 + 벤치마크 확장

## Topic B: ssPMT layout 공개의 leakage 분석과 보수

### 문제 설정
KKLS ssPMT를 단독 사용할 때, 성능 최적화로 "R의 점유 decode position(layout/bitset)을
S에게 공개"하면 ss-equality를 L·H개가 아닌 n_y개 인스턴스로 압축할 수 있다. 이 공개가
안전한가?

### 결론: 조건부로 불안전 (누출이 λ = n_y/H에 의해 지배됨)

**공격 전제**: S는 indicator OKVS를 F(x)로 인코딩하므로 자기 아이템의 F(x)를 알고,
따라서 slot 좌표 τ(F(x)) (band start H₁(F(x))에서 결정)를 예측할 수 있다.

**공격 (slot-column 다중도 LR 테스트)**:
- 관찰량: 열 τ의 점유수 m_τ = |{j : τ(F(y_j)) = τ}|.
- 분포: 비멤버 x → m_{τ(F(x))} ~ Bin(n_y, 1/H) ≈ Poi(λ); 멤버 x → 1 + Poi(λ).
- 우도비 LR(m) = m/λ → **Bayes 최적 판정: m ≥ λ+1 이면 멤버로 판정** (LR ≥ 1 ⟺ m ≥ λ).
- 후보당 구별 우위 = TV(Poi(λ), 1+Poi(λ)) = **P(mode) ≈ 1/√(2πλ)**
  (단봉 분포의 한 칸 shift에 대한 telescoping 항등식: TV = 최빈값 질량).
  λ=128이면 ≈3.5%, 성공확률 ≈ ½ + 1.8%.
- **질적으로 위험한 부분: m=0 확정 배제 인증서.** LR(0)=0이므로 m=0이면 비멤버 확정.
  발생확률 e^{−λ}: well-filled(λ=128)에선 사실상 0이지만, **성긴 layout(λ~1)이면
  비멤버의 ~37%가 확정 배제** — 배치 파라미터에 따라 질적으로 다른 누출.
- 집계 공격(Σm로 |X∩Y| 추정)은 교집합 크기가 어차피 functionality 출력이라 새 누출 아님.

**더 나쁜 변형 (중요)**: 다른 구조(예: circuit-PSI의 cuckoo bin)와의 **대응표
σ: bin ↦ (μ,τ)를 공개하면**, S가 (bin(F(x)), τ(F(x))) 쌍의 일관성을 확인하는
consistency oracle로 **오차 ~1/H의 사실상 완전한 per-item 멤버십 테스트**가 됨.
Dummy padding으로 못 막음 (real 아이템의 쌍은 여전히 일관). S가 예측 가능한 두 좌표의
짝을 공개하는 순간 죽는다.

### 보수 (권장 순서)
1. **공개하지 않기 (기본값)**: 정해진(deterministic) full layout L·H 전체에 대해
   ss-equality 수행. 비용 = 1/fill-rate ≈ 1.1–1.4× (상수 배). 또는 동치로 R이 dummy
   질의로 전 slot을 채움. S의 view에서 Y-의존 정보가 소멸 → 시뮬레이션 자명.
2. 압축이 꼭 필요하면: leakage function을 명시하고 (위 다중도 분포 + TV 정량화) 그
   leakage 하의 보안으로 증명. sparse 배치 금지 조건(λ 하한) 명시.
3. 부분 완화(예: layer별 비밀 rotation으로 열 좌표 은닉)는 rotation-invariant한 gap
   구조가 남아 pair-correlation 공격 가능 → 깨끗하지 않음, 비권장.

### 후속 논문이 할 일
1. 위 공격의 정식화 (leakage function, 최적성: Neyman–Pearson으로 LR 판정 최적 증명)
2. TV = P(mode) 항등식 정리로 정리(단봉 shift 일반), λ-스케일링 명시
3. Full-layout(또는 dummy-fill) 버전의 시뮬레이션 증명
4. 비용 실측: n_y vs L·H 인스턴스의 ss-equality 비용 차이
5. (검증 필요) consistency-oracle 변형의 서술 — 원 세션 유도라 독립 재검토 권장

## 하지 말 것 / 경계

- σ_min tail bound, real-field OKVS, min-norm encoding 등은 **Real-PJC 논문 소관** —
  이 후속 연구 범위 아님.
- Real-PJC 쪽에서도 같은 leakage 논의가 "unified layout" 논증으로 들어갈 예정이지만,
  거기서는 alignment 프로토콜의 설계 정당화 용도이고, 여기서는 KKLS ssPMT 자체의
  분석·보수가 목적. 서술 중복은 괜찮으나 결론이 어긋나지 않게 유지할 것.

## 시작 프롬프트 예시

> docs/kkls-followup-brief.md를 읽고 KKLS 후속 note 작업을 시작해줘. 우선 (1) design
> note에서 dense two-row 관련 부분과 git 히스토리(2d828a1^)의 삭제된 섹션·appendix를
> 회수해 새 tex 골격을 만들고, (2) Topic B의 leakage 분석을 정식 서술로 옮겨줘.
