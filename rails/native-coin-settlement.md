# Native-Coin Settlement — Nano (XNO) under x402

> A fee-free, instantly-final native-coin settlement rail for agent-to-API payments. Where stablecoin rails carry an issuer (USDC/USDT/DAI/EURC) and a gas fee on every transfer, a native-coin rail settles the protocol's own currency with zero transaction fees and no issuer to freeze or gate a balance. This page covers how a merchant accepts it today over the x402 framing.

## What this is

Native-coin settlement means the buyer pays in the settlement network's own currency — for Nano, **XNO** — rather than in a stablecoin issued by an external party. The transfer is a single block on a block-lattice DAG; there is no gas, no fee market, and no issuer-mint/burn leg. From an agent-buyer perspective the shape is unchanged from stablecoin x402: the server returns an HTTP 402 with a payment challenge, the agent signs and broadcasts a transfer, and the retried request carries the proof.

Mechanically:

1. The agent sends the request and receives `HTTP 402 Payment Required` with a challenge whose `accepts` names the settlement asset (e.g. `scheme: exact`, `network: nano:mainnet`, `asset: XNO`, `amount`, `payTo`).
2. The agent's wallet builds a signed Nano send block to `payTo` for the exact amount.
3. The block is broadcast and confirmed on the Nano network.
4. The merchant detects the on-chain balance change and matches it to the order via an idempotency key carried in the quote.
5. The merchant releases the resource.

## Why a merchant would add it

- **Zero transaction fees.** Nano is fee-less by protocol design — there is no per-transfer cost for either party, so sub-cent unit economics are not capped by a gas floor.
- **Instant, final settlement.** Nano uses a block-lattice DAG with ~0.3 s block times and no reorg probability, so a merchant does not wait out a probabilistic finality window (compare Base's ~12 min hard finality or Bitcoin's 6 confirmations).
- **No issuer risk.** The asset is not an IOU from an issuer that can freeze or blacklist a balance; the merchant holds the protocol's own currency.
- **True micropayments.** Because there is no fee to amortize, per-call pricing down to fractions of a cent is practical without batching or channel topology.

## What it does not change

- The dispute model is the same as other crypto rails: once the block is on-chain it is final; refunds are merchant-initiated reverse transfers. Encode refund policy in the receipt (see `/merchant-playbooks/refunds-and-disputes-for-agents.md`).
- Idempotency keys, cost ceilings, and quota management apply exactly as they do to stablecoin x402 (see `/use-cases/m2m-and-api.md`).

## Live example

A pay-per-call API settled in Nano under the x402 framing is Vend's API merchant. Its manifest advertises the `accepts` payload with `network: nano:mainnet`, `asset: XNO`, a 0.0001 XNO price and a treasury `payTo`, and returns HTTP 200 at:

- <https://extract.paypercall.dev/.well-known/x402>

The same change is documented in the cross-rail comparison at `/comparison/rails-comparison.md`.

## References

- Nano docs, "What is nano?" — zero fees, ultrafast, "cannot be reversed" — <https://docs.nano.org/what-is-nano/overview/>
- Nano protocol design — <https://docs.nano.org/protocol-design/introduction/>
- x402 specification — <https://www.x402.org/>
- Vend x402 manifest (live, HTTP 200) — <https://extract.paypercall.dev/.well-known/x402>
- Cross-rail comparison — [/comparison/rails-comparison.md](../comparison/rails-comparison.md)
