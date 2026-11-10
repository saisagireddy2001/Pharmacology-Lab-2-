"""Reproduce analysis plots and summary tables for the labelled teaching CSV data."""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = {'time_s', 'twitch_force_arbitrary_units', 'target_twitch_percent_baseline', 'stage', 'stimulus', 'dataset_type'}
SUXA_UM = [0.2, 0.4, 1.5, 6, 12, 25, 50, 100]


def response_percent(twitch: float, baseline: float) -> float:
    if baseline <= 0:
        raise ValueError('Baseline force must be positive')
    return 100.0 * twitch / baseline


def inhibition_percent(twitch: float, baseline: float) -> float:
    return 100.0 - response_percent(twitch, baseline)


def log_molar(micromolar: float) -> float:
    if micromolar <= 0:
        raise ValueError('Concentration must be positive')
    return float(np.log10(micromolar * 1e-6))


def analyze_csv(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(path)
    missing = REQUIRED.difference(df.columns)
    if missing:
        raise ValueError(f'{path.name}: missing columns {sorted(missing)}')
    if not df['dataset_type'].eq('SYNTHETIC_SIMULATION').all():
        raise ValueError(f'{path.name}: unexpected provenance; use validated parser for other data')
    # group contiguous stages; repeated 'Wash' labels must remain separate events
    stage_num = df['stage'].ne(df['stage'].shift()).cumsum()
    summary = df.groupby(stage_num, sort=True).agg(
        stage=('stage','first'),
        start_s=('time_s','min'),
        end_s=('time_s','max'),
        response_pct=('target_twitch_percent_baseline','first'),
        n_samples=('time_s','size')
    ).reset_index(drop=True)
    summary['inhibition_pct'] = 100 - summary['response_pct']
    summary.insert(0,'experiment',path.stem)
    return df, summary


def make_trace_plot(df: pd.DataFrame, summary: pd.DataFrame, target: Path) -> None:
    fig, ax = plt.subplots(figsize=(13, 4))
    ax.plot(df.time_s, df.twitch_force_arbitrary_units, lw=.65)
    for _, s in summary.iterrows():
        ax.axvline(s.start_s, color='gray', lw=.6, alpha=.4)
    ax.set(title=f'{target.stem.replace("_", " ")} — simulated trace', xlabel='Time (s)', ylabel='Force (arbitrary units)')
    fig.tight_layout()
    fig.savefig(target, dpi=140)
    plt.close(fig)


def make_suxa_table(summary: pd.DataFrame) -> pd.DataFrame:
    measured = summary[summary.stage.str.startswith('Suxamethonium ')].copy()
    if len(measured) != 8:
        raise ValueError('Expected eight suxamethonium dosing stages')
    measured['concentration_uM'] = SUXA_UM
    measured['concentration_M'] = measured.concentration_uM * 1e-6
    measured['log10_concentration_M'] = measured.concentration_uM.map(log_molar)
    return measured[['concentration_uM','concentration_M','log10_concentration_M','response_pct','inhibition_pct']]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--input-dir', type=Path, default=ROOT / 'data/simulated_csv')
    p.add_argument('--output-dir', type=Path, default=ROOT / 'outputs')
    p.add_argument('--no-trace-plots', action='store_true')
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_summaries=[]
    suxa=None
    files=sorted(args.input_dir.glob('*.csv'))
    if not files:
        raise SystemExit('No CSV trace files found')
    for f in files:
        df, stages=analyze_csv(f)
        all_summaries.append(stages)
        if not args.no_trace_plots:
            make_trace_plot(df,stages,args.output_dir/(f.stem+'.png'))
        if f.stem == 'Part_III_Suxamethonium_Eight_Dose_Curve':
            suxa=make_suxa_table(stages)
    pd.concat(all_summaries,ignore_index=True).to_csv(args.output_dir/'stage_summary.csv',index=False)
    if suxa is not None:
        suxa.to_csv(args.output_dir/'suxamethonium_response.csv',index=False)
        fig,ax=plt.subplots(figsize=(7.5,4.5))
        ax.plot(suxa.log10_concentration_M,suxa.inhibition_pct,'o-',lw=2)
        ax.set(xlabel='log10 [suxamethonium] (M)', ylabel='Twitch inhibition (%)',title='Suxamethonium concentration–response (teaching model)')
        ax.set_ylim(-5,105)
        fig.tight_layout()
        fig.savefig(args.output_dir/'suxamethonium_dose_response.png',dpi=160)
        plt.close(fig)
    print(f'Analyzed {len(files)} CSV traces. Output: {args.output_dir}')

if __name__ == '__main__':
    main()
