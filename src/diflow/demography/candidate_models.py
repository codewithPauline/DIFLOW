"""Candidate two-population demographic models for DIFLOW.

All functions return dadi-compatible spectra while DIFLOW keeps forward-time
source-to-recipient biological labels for migration parameters.
"""

from __future__ import annotations


def no_migration_model(dadi):
    def model(params, ns, pts):
        nu_a, nu_b, split_time = params
        return dadi.Demographics2D.no_mig(
            (nu_a, nu_b, split_time),
            ns,
            pts,
        )
    return model


def symmetric_migration_model(dadi):
    def model(params, ns, pts):
        nu_a, nu_b, split_time, migration = params
        return dadi.Demographics2D.split_mig(
            (nu_a, nu_b, split_time, migration),
            ns,
            pts,
        )
    return model


def asymmetric_migration_model(dadi):
    def model(params, ns, pts):
        nu_a, nu_b, split_time, m_a_to_b, m_b_to_a = params
        return dadi.Demographics2D.split_asym_mig(
            (
                nu_a,
                nu_b,
                split_time,
                m_b_to_a,
                m_a_to_b,
            ),
            ns,
            pts,
        )
    return model


def secondary_contact_asymmetric_model(dadi):
    def model(params, ns, pts):
        nu_a, nu_b, isolation_time, contact_time, m_a_to_b, m_b_to_a = params
        return dadi.Demographics2D.sec_contact_asym_mig(
            (
                nu_a,
                nu_b,
                m_b_to_a,
                m_a_to_b,
                isolation_time,
                contact_time,
            ),
            ns,
            pts,
        )
    return model
