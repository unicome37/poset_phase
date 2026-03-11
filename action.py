from __future__ import annotations

from observables import neutral_penalty
from observables_geo import geometric_penalty


def action_value(log_extensions: float, penalty: float, beta: float, gamma: float) -> float:
    return -beta * log_extensions + gamma * penalty


def score_neutral(poset, log_extensions: float, beta: float, gamma: float) -> float:
    return action_value(
        log_extensions=log_extensions,
        penalty=neutral_penalty(poset),
        beta=beta,
        gamma=gamma,
    )


def score_geometric(poset, log_extensions: float, beta: float, gamma: float) -> float:
    return action_value(
        log_extensions=log_extensions,
        penalty=geometric_penalty(poset),
        beta=beta,
        gamma=gamma,
    )


def score_mixed(poset, log_extensions: float, beta: float, gamma: float) -> float:
    return action_value(
        log_extensions=log_extensions,
        penalty=neutral_penalty(poset) + geometric_penalty(poset),
        beta=beta,
        gamma=gamma,
    )


def get_action_penalty(poset, action_mode: str) -> float:
    if action_mode == "A1":
        return neutral_penalty(poset)
    if action_mode == "A2":
        return neutral_penalty(poset) + geometric_penalty(poset)
    if action_mode == "A3":
        return geometric_penalty(poset)
    raise ValueError(f"Unsupported action mode: {action_mode}")
