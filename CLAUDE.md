# CLAUDE.md

이 레포는 KKLS(base RPMT, EPRINT:KwaKwoLeeSon26)의 구현(`rlwe-okvs/`, `exec/`,
`experiments/`)을 담는다. (후속 note와 CPSI 코드는 2026-08-12에 제거됨 —
git 히스토리에는 남아 있다.)

## 사용자 선호 — 채팅 응답에서의 수식 표기

- 채팅/마크다운 응답에서 수식을 LaTeX 문법(`$...$`, `\mu`, `\sqrt{}` 등)으로 쓰지
  말 것. 대신 유니코드로 쓸 것: μ, τ, λ = n_y/H, e^−λ, 1/√(2πλ), ⌊λ⌋, Σ, ⊕, ≈, ≤ 등.
- 이 규칙은 채팅 응답에만 적용된다. `.tex` 파일 등 문서 작성 시에는 LaTeX를 쓴다.
