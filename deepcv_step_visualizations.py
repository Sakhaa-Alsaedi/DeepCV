#!/usr/bin/env python3
"""
DeepCV Step-by-Step Visualizations Module

This module provides comprehensive visualizations for each step of the DeepCV pipeline:
- Step 1: VCF Parsing and Quality Control
- Step 2: SNP Annotation Analysis
- Step 3: Risk Assessment Visualizations
- Step 4: Network Analysis (already implemented)

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class DeepCVStepVisualizations:
    """
    Comprehensive visualization class for DeepCV pipeline steps
    """
    
    def __init__(self, output_dir="results"):
        """
        Initialize the visualization class
        
        Args:
            output_dir (str): Directory to save visualizations
        """
        self.output_dir = Path(output_dir)
        self.viz_dir = self.output_dir / "step_visualizations"
        self.viz_dir.mkdir(parents=True, exist_ok=True)
        
        # Set plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
    def create_step1_visualizations(self, vcf_data_file, vcf_basename=""):
        """
        Create comprehensive visualizations for Step 1: VCF Parsing
        
        Args:
            vcf_data_file (str): Path to processed variants CSV file
            vcf_basename (str): Base name for output files
            
        Returns:
            dict: Dictionary of created plot paths
        """
        print(f"\n📊 Creating Step 1 Visualizations: VCF Parsing & Quality Control")
        
        try:
            # Load VCF data
            df = pd.read_csv(vcf_data_file)
            print(f"   📁 Loaded {len(df)} variants from {vcf_data_file}")
            
            plots_created = {}
            
            # 1. Quality Score Distributions
            quality_plot = self._create_quality_distributions(df, vcf_basename)
            if quality_plot:
                plots_created['quality_distributions'] = quality_plot
            
            # 2. Manhattan-like Plot
            manhattan_plot = self._create_manhattan_plot(df, vcf_basename)
            if manhattan_plot:
                plots_created['manhattan_plot'] = manhattan_plot
            
            # 3. Variant Type Distribution
            variant_type_plot = self._create_variant_type_distribution(df, vcf_basename)
            if variant_type_plot:
                plots_created['variant_types'] = variant_type_plot
            
            # 4. Chromosome Distribution
            chr_dist_plot = self._create_chromosome_distribution(df, vcf_basename)
            if chr_dist_plot:
                plots_created['chromosome_distribution'] = chr_dist_plot
            
            print(f"   Created {len(plots_created)} Step 1 visualizations")
            return plots_created
            
        except Exception as e:
            print(f"   Error creating Step 1 visualizations: {e}")
            return {}
    
    def _create_quality_distributions(self, df, vcf_basename):
        """Create quality score distribution plots"""
        try:
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            fig.suptitle('VCF Quality Score Distributions', fontsize=16, fontweight='bold')
            
            # VAF Distribution
            if 'VAF' in df.columns:
                axes[0, 0].hist(df['VAF'].dropna(), bins=50, alpha=0.7, color='skyblue', edgecolor='black')
                axes[0, 0].set_xlabel('Variant Allele Frequency (VAF)')
                axes[0, 0].set_ylabel('Frequency')
                axes[0, 0].set_title('VAF Distribution')
                axes[0, 0].grid(True, alpha=0.3)
            
            # GQ Distribution
            if 'GQ' in df.columns:
                axes[0, 1].hist(df['GQ'].dropna(), bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
                axes[0, 1].set_xlabel('Genotype Quality (GQ)')
                axes[0, 1].set_ylabel('Frequency')
                axes[0, 1].set_title('GQ Distribution')
                axes[0, 1].grid(True, alpha=0.3)
            
            # QUAL Distribution
            if 'QUAL' in df.columns:
                axes[1, 0].hist(df['QUAL'].dropna(), bins=50, alpha=0.7, color='salmon', edgecolor='black')
                axes[1, 0].set_xlabel('Quality Score (QUAL)')
                axes[1, 0].set_ylabel('Frequency')
                axes[1, 0].set_title('QUAL Distribution')
                axes[1, 0].grid(True, alpha=0.3)
            
            # DP Distribution (if available)
            if 'DP' in df.columns:
                axes[1, 1].hist(df['DP'].dropna(), bins=50, alpha=0.7, color='gold', edgecolor='black')
                axes[1, 1].set_xlabel('Read Depth (DP)')
                axes[1, 1].set_ylabel('Frequency')
                axes[1, 1].set_title('DP Distribution')
                axes[1, 1].grid(True, alpha=0.3)
            else:
                # Alternative: Position distribution
                if 'POS' in df.columns:
                    axes[1, 1].hist(df['POS'].dropna(), bins=50, alpha=0.7, color='purple', edgecolor='black')
                    axes[1, 1].set_xlabel('Position')
                    axes[1, 1].set_ylabel('Frequency')
                    axes[1, 1].set_title('Position Distribution')
                    axes[1, 1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            output_path = self.viz_dir / f"step1_quality_distributions{vcf_basename}.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"   Quality distributions saved: {output_path}")
            return str(output_path)
            
        except Exception as e:
            print(f"   Error creating quality distributions: {e}")
            return None
    
    def _create_manhattan_plot(self, df, vcf_basename):
        """Create Manhattan-like plot for variant positions"""
        try:
            fig, ax = plt.subplots(1, 1, figsize=(16, 8))
            
            # Prepare data for Manhattan plot
            if 'CHROM' in df.columns and 'POS' in df.columns:
                # Convert chromosome to numeric for plotting
                df_plot = df.copy()
                
                # Handle chromosome naming (remove 'chr' prefix if present)
                df_plot['CHROM_clean'] = df_plot['CHROM'].astype(str).str.replace('chr', '')
                
                # Convert to numeric, handling X, Y, MT
                chrom_map = {'X': 23, 'Y': 24, 'MT': 25, 'M': 25}
                df_plot['CHROM_num'] = df_plot['CHROM_clean'].replace(chrom_map)
                df_plot['CHROM_num'] = pd.to_numeric(df_plot['CHROM_num'], errors='coerce')
                
                # Remove invalid chromosomes
                df_plot = df_plot.dropna(subset=['CHROM_num'])
                
                # Create cumulative positions for x-axis
                df_plot = df_plot.sort_values(['CHROM_num', 'POS'])
                
                # Calculate cumulative positions
                chrom_lengths = df_plot.groupby('CHROM_num')['POS'].max()
                cumulative_pos = 0
                df_plot['cumulative_pos'] = 0
                
                for chrom in sorted(df_plot['CHROM_num'].unique()):
                    mask = df_plot['CHROM_num'] == chrom
                    df_plot.loc[mask, 'cumulative_pos'] = df_plot.loc[mask, 'POS'] + cumulative_pos
                    cumulative_pos += chrom_lengths[chrom]
                
                # Use VAF for y-axis if available, otherwise use a constant
                if 'VAF' in df_plot.columns:
                    y_values = df_plot['VAF']
                    y_label = 'Variant Allele Frequency (VAF)'
                else:
                    y_values = np.ones(len(df_plot))  # Constant height
                    y_label = 'Variant Presence'
                
                # Create color map for chromosomes
                colors = plt.cm.Set3(np.linspace(0, 1, len(df_plot['CHROM_num'].unique())))
                
                for i, chrom in enumerate(sorted(df_plot['CHROM_num'].unique())):
                    mask = df_plot['CHROM_num'] == chrom
                    ax.scatter(df_plot.loc[mask, 'cumulative_pos'], 
                             y_values[mask], 
                             c=[colors[i]], 
                             alpha=0.6, 
                             s=20, 
                             label=f'Chr {int(chrom) if chrom <= 22 else {23: "X", 24: "Y", 25: "MT"}[chrom]}')
                
                ax.set_xlabel('Genomic Position (bp)')
                ax.set_ylabel(y_label)
                ax.set_title('Manhattan Plot: Variant Distribution Across Chromosomes')
                ax.grid(True, alpha=0.3)
                
                # Add chromosome boundaries
                chrom_centers = []
                cumulative_pos = 0
                for chrom in sorted(df_plot['CHROM_num'].unique()):
                    chrom_length = chrom_lengths[chrom]
                    chrom_center = cumulative_pos + chrom_length / 2
                    chrom_centers.append(chrom_center)
                    cumulative_pos += chrom_length
                    
                    # Add vertical line at chromosome boundary
                    if chrom < max(df_plot['CHROM_num']):
                        ax.axvline(x=cumulative_pos, color='gray', linestyle='--', alpha=0.5)
                
                # Set x-axis labels to chromosome names
                ax.set_xticks(chrom_centers)
                ax.set_xticklabels([f'{int(c) if c <= 22 else {23: "X", 24: "Y", 25: "MT"}[c]}' 
                                   for c in sorted(df_plot['CHROM_num'].unique())])
                
                plt.tight_layout()
                
                output_path = self.viz_dir / f"step1_manhattan_plot{vcf_basename}.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   Manhattan plot saved: {output_path}")
                return str(output_path)
            
        except Exception as e:
            print(f"   Error creating Manhattan plot: {e}")
            return None
    
    def _create_variant_type_distribution(self, df, vcf_basename):
        """Create variant type distribution plot"""
        try:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            
            # Determine variant types based on REF and ALT
            if 'REF' in df.columns and 'ALT' in df.columns:
                df_clean = df.dropna(subset=['REF', 'ALT'])
                
                variant_types = []
                for _, row in df_clean.iterrows():
                    ref = str(row['REF'])
                    alt = str(row['ALT'])
                    
                    if len(ref) == 1 and len(alt) == 1:
                        variant_types.append('SNV')
                    elif len(ref) > len(alt):
                        variant_types.append('Deletion')
                    elif len(ref) < len(alt):
                        variant_types.append('Insertion')
                    else:
                        variant_types.append('Complex')
                
                # Variant type counts
                type_counts = pd.Series(variant_types).value_counts()
                
                # Bar plot
                ax1.bar(type_counts.index, type_counts.values, alpha=0.7, color=['skyblue', 'lightgreen', 'salmon', 'gold'])
                ax1.set_xlabel('Variant Type')
                ax1.set_ylabel('Count')
                ax1.set_title('Variant Type Distribution')
                ax1.grid(True, alpha=0.3)
                
                # Pie chart
                ax2.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%', startangle=90)
                ax2.set_title('Variant Type Proportions')
            
            plt.tight_layout()
            
            output_path = self.viz_dir / f"step1_variant_types{vcf_basename}.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"   Variant types plot saved: {output_path}")
            return str(output_path)
            
        except Exception as e:
            print(f"   Error creating variant types plot: {e}")
            return None
    
    def _create_chromosome_distribution(self, df, vcf_basename):
        """Create chromosome distribution plot"""
        try:
            if 'CHROM' in df.columns:
                fig, ax = plt.subplots(1, 1, figsize=(14, 8))
                
                # Count variants per chromosome
                chrom_counts = df['CHROM'].value_counts().sort_index()
                
                # Create bar plot
                bars = ax.bar(range(len(chrom_counts)), chrom_counts.values, alpha=0.7)
                
                # Color bars alternately for better visualization
                for i, bar in enumerate(bars):
                    if i % 2 == 0:
                        bar.set_color('skyblue')
                    else:
                        bar.set_color('lightcoral')
                
                ax.set_xlabel('Chromosome')
                ax.set_ylabel('Number of Variants')
                ax.set_title('Variant Distribution Across Chromosomes')
                ax.set_xticks(range(len(chrom_counts)))
                ax.set_xticklabels(chrom_counts.index, rotation=45)
                ax.grid(True, alpha=0.3)
                
                # Add value labels on bars
                for i, v in enumerate(chrom_counts.values):
                    ax.text(i, v + max(chrom_counts.values) * 0.01, str(v), 
                           ha='center', va='bottom', fontsize=9)
                
                plt.tight_layout()
                
                output_path = self.viz_dir / f"step1_chromosome_distribution{vcf_basename}.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   Chromosome distribution saved: {output_path}")
                return str(output_path)
                
        except Exception as e:
            print(f"   Error creating chromosome distribution: {e}")
            return None
    
    def create_step2_visualizations(self, annotated_file, original_file, vcf_basename=""):
        """
        Create visualizations for Step 2: SNP Annotation
        
        Args:
            annotated_file (str): Path to annotated variants file
            original_file (str): Path to original variants file
            vcf_basename (str): Base name for output files
            
        Returns:
            dict: Dictionary of created plot paths and annotation summary
        """
        print(f"\n Creating Step 2 Visualizations: SNP Annotation Analysis")
        
        try:
            # Load data
            df_annotated = pd.read_csv(annotated_file)
            df_original = pd.read_csv(original_file)
            
            print(f"   📁 Original: {len(df_original)} variants, Annotated: {len(df_annotated)} variants")
            
            plots_created = {}
            
            # 1. Annotation Features Summary
            annotation_summary = self._create_annotation_summary(df_annotated, df_original, vcf_basename)
            if annotation_summary:
                plots_created['annotation_summary'] = annotation_summary
            
            # 2. Gene Impact Distribution
            impact_plot = self._create_impact_distribution(df_annotated, vcf_basename)
            if impact_plot:
                plots_created['impact_distribution'] = impact_plot
            
            # 3. Consequence Type Distribution
            consequence_plot = self._create_consequence_distribution(df_annotated, vcf_basename)
            if consequence_plot:
                plots_created['consequence_distribution'] = consequence_plot
            
            # 4. Annotation Coverage Analysis
            coverage_plot = self._create_annotation_coverage(df_annotated, vcf_basename)
            if coverage_plot:
                plots_created['annotation_coverage'] = coverage_plot
            
            print(f"   Created {len(plots_created)} Step 2 visualizations")
            return plots_created
            
        except Exception as e:
            print(f"   Error creating Step 2 visualizations: {e}")
            return {}
    
    def _create_annotation_summary(self, df_annotated, df_original, vcf_basename):
        """Create annotation summary comparison"""
        try:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
            fig.suptitle('SNP Annotation Summary', fontsize=16, fontweight='bold')
            
            # 1. Column comparison
            original_cols = set(df_original.columns)
            annotated_cols = set(df_annotated.columns)
            new_cols = annotated_cols - original_cols
            
            col_data = {
                'Original Columns': len(original_cols),
                'New Annotations': len(new_cols),
                'Total Columns': len(annotated_cols)
            }
            
            ax1.bar(col_data.keys(), col_data.values(), color=['lightblue', 'lightgreen', 'salmon'])
            ax1.set_ylabel('Number of Columns')
            ax1.set_title('Annotation Features Added')
            ax1.grid(True, alpha=0.3)
            
            # Add value labels
            for i, (k, v) in enumerate(col_data.items()):
                ax1.text(i, v + max(col_data.values()) * 0.01, str(v), 
                        ha='center', va='bottom', fontweight='bold')
            
            # 2. New annotation features list
            if new_cols:
                new_cols_list = sorted(list(new_cols))
                y_pos = np.arange(len(new_cols_list))
                
                ax2.barh(y_pos, [1] * len(new_cols_list), color='lightgreen', alpha=0.7)
                ax2.set_yticks(y_pos)
                ax2.set_yticklabels(new_cols_list, fontsize=8)
                ax2.set_xlabel('Feature Added')
                ax2.set_title(f'New Annotation Features ({len(new_cols_list)})')
                ax2.grid(True, alpha=0.3)
            
            # 3. Data completeness comparison
            original_completeness = (df_original.notna().sum() / len(df_original) * 100).mean()
            annotated_completeness = (df_annotated.notna().sum() / len(df_annotated) * 100).mean()
            
            completeness_data = {
                'Original Data': original_completeness,
                'Annotated Data': annotated_completeness
            }
            
            ax3.bar(completeness_data.keys(), completeness_data.values(), 
                   color=['lightcoral', 'lightgreen'])
            ax3.set_ylabel('Data Completeness (%)')
            ax3.set_title('Data Completeness Comparison')
            ax3.set_ylim(0, 100)
            ax3.grid(True, alpha=0.3)
            
            # Add value labels
            for i, (k, v) in enumerate(completeness_data.items()):
                ax3.text(i, v + 2, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')
            
            # 4. Annotation success rate
            if 'Gene' in df_annotated.columns:
                annotated_genes = df_annotated['Gene'].notna().sum()
                total_variants = len(df_annotated)
                success_rate = (annotated_genes / total_variants) * 100
                
                ax4.pie([success_rate, 100 - success_rate], 
                       labels=[f'Annotated ({annotated_genes})', f'Not Annotated ({total_variants - annotated_genes})'],
                       autopct='%1.1f%%', startangle=90, colors=['lightgreen', 'lightcoral'])
                ax4.set_title('Gene Annotation Success Rate')
            
            plt.tight_layout()
            
            output_path = self.viz_dir / f"step2_annotation_summary{vcf_basename}.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"   Annotation summary saved: {output_path}")
            
            # Save new features list as text file
            if new_cols:
                features_file = self.viz_dir / f"step2_new_features{vcf_basename}.txt"
                with open(features_file, 'w') as f:
                    f.write("New Annotation Features Added in Step 2:\n")
                    f.write("=" * 50 + "\n\n")
                    for i, feature in enumerate(sorted(new_cols), 1):
                        f.write(f"{i:2d}. {feature}\n")
                    f.write(f"\nTotal new features: {len(new_cols)}\n")
                
                print(f"   New features list saved: {features_file}")
            
            return str(output_path)
            
        except Exception as e:
            print(f"   Error creating annotation summary: {e}")
            return None
    
    def _create_impact_distribution(self, df, vcf_basename):
        """Create variant impact distribution plot"""
        try:
            # Look for impact-related columns
            impact_cols = [col for col in df.columns if 'impact' in col.lower() or 'severity' in col.lower()]
            
            if impact_cols:
                fig, axes = plt.subplots(1, len(impact_cols), figsize=(6 * len(impact_cols), 6))
                if len(impact_cols) == 1:
                    axes = [axes]
                
                for i, col in enumerate(impact_cols):
                    impact_counts = df[col].value_counts()
                    
                    # Create pie chart
                    axes[i].pie(impact_counts.values, labels=impact_counts.index, 
                               autopct='%1.1f%%', startangle=90)
                    axes[i].set_title(f'{col} Distribution')
                
                plt.tight_layout()
                
                output_path = self.viz_dir / f"step2_impact_distribution{vcf_basename}.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   🎯 Impact distribution saved: {output_path}")
                return str(output_path)
                
        except Exception as e:
            print(f"   Error creating impact distribution: {e}")
            return None
    
    def _create_consequence_distribution(self, df, vcf_basename):
        """Create consequence type distribution plot"""
        try:
            # Look for consequence-related columns
            consequence_cols = [col for col in df.columns if 'consequence' in col.lower() or 'effect' in col.lower()]
            
            if consequence_cols:
                fig, ax = plt.subplots(1, 1, figsize=(12, 8))
                
                col = consequence_cols[0]  # Use first consequence column
                consequence_counts = df[col].value_counts().head(15)  # Top 15
                
                bars = ax.bar(range(len(consequence_counts)), consequence_counts.values, 
                             color=plt.cm.Set3(np.linspace(0, 1, len(consequence_counts))))
                
                ax.set_xlabel('Consequence Type')
                ax.set_ylabel('Count')
                ax.set_title(f'Top 15 {col} Types')
                ax.set_xticks(range(len(consequence_counts)))
                ax.set_xticklabels(consequence_counts.index, rotation=45, ha='right')
                ax.grid(True, alpha=0.3)
                
                # Add value labels
                for i, v in enumerate(consequence_counts.values):
                    ax.text(i, v + max(consequence_counts.values) * 0.01, str(v), 
                           ha='center', va='bottom', fontsize=8)
                
                plt.tight_layout()
                
                output_path = self.viz_dir / f"step2_consequence_distribution{vcf_basename}.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   Consequence distribution saved: {output_path}")
                return str(output_path)
                
        except Exception as e:
            print(f"   Error creating consequence distribution: {e}")
            return None
    
    def _create_annotation_coverage(self, df, vcf_basename):
        """Create annotation coverage heatmap"""
        try:
            # Calculate missing data percentage for each column
            missing_data = (df.isnull().sum() / len(df) * 100).sort_values(ascending=True)
            
            # Select annotation columns (exclude basic VCF columns)
            basic_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
            annotation_cols = [col for col in missing_data.index if col not in basic_cols]
            
            if annotation_cols:
                fig, ax = plt.subplots(1, 1, figsize=(12, max(8, len(annotation_cols) * 0.3)))
                
                # Create horizontal bar plot
                y_pos = np.arange(len(annotation_cols))
                missing_values = missing_data[annotation_cols].values
                
                bars = ax.barh(y_pos, 100 - missing_values, color='lightgreen', alpha=0.7, label='Available')
                ax.barh(y_pos, missing_values, left=100 - missing_values, color='lightcoral', alpha=0.7, label='Missing')
                
                ax.set_yticks(y_pos)
                ax.set_yticklabels(annotation_cols, fontsize=8)
                ax.set_xlabel('Percentage')
                ax.set_title('Annotation Coverage by Feature')
                ax.legend()
                ax.grid(True, alpha=0.3)
                
                plt.tight_layout()
                
                output_path = self.viz_dir / f"step2_annotation_coverage{vcf_basename}.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   Annotation coverage saved: {output_path}")
                return str(output_path)
                
        except Exception as e:
            print(f"   Error creating annotation coverage: {e}")
            return None
    
    def create_step3_visualizations(self, risk_assessment_file, vcf_basename=""):
        """
        Create visualizations for Step 3: Risk Assessment
        
        Args:
            risk_assessment_file (str): Path to risk assessment file
            vcf_basename (str): Base name for output files
            
        Returns:
            dict: Dictionary of created plot paths
        """
        print(f"\n Creating Step 3 Visualizations: Risk Assessment Analysis")
        
        try:
            # Load risk assessment data
            df = pd.read_csv(risk_assessment_file)
            print(f"   Loaded {len(df)} risk associations from {risk_assessment_file}")
            
            plots_created = {}
            
            # 1. Risk Score Distribution
            risk_dist_plot = self._create_risk_score_distribution(df, vcf_basename)
            if risk_dist_plot:
                plots_created['risk_distribution'] = risk_dist_plot
            
            # 2. Top Genes Analysis
            top_genes_plot = self._create_top_genes_analysis(df, vcf_basename)
            if top_genes_plot:
                plots_created['top_genes'] = top_genes_plot
            
            # 3. Disease Categories Analysis
            disease_plot = self._create_disease_analysis(df, vcf_basename)
            if disease_plot:
                plots_created['disease_analysis'] = disease_plot
            
            # 4. Gene-Disease Heatmap
            heatmap_plot = self._create_gene_disease_heatmap(df, vcf_basename)
            if heatmap_plot:
                plots_created['gene_disease_heatmap'] = heatmap_plot
            
            print(f"   Created {len(plots_created)} Step 3 visualizations")
            return plots_created
            
        except Exception as e:
            print(f"   Error creating Step 3 visualizations: {e}")
            return {}
    
    def _create_risk_score_distribution(self, df, vcf_basename):
        """Create risk score distribution analysis"""
        try:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
            fig.suptitle('Risk Score Distribution Analysis', fontsize=16, fontweight='bold')
            
            # 1. Overall risk score distribution
            if 'score' in df.columns:
                ax1.hist(df['score'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
                ax1.set_xlabel('Risk Score')
                ax1.set_ylabel('Frequency')
                ax1.set_title('Risk Score Distribution')
                ax1.grid(True, alpha=0.3)
                
                # Add statistics
                mean_score = df['score'].mean()
                median_score = df['score'].median()
                ax1.axvline(mean_score, color='red', linestyle='--', label=f'Mean: {mean_score:.3f}')
                ax1.axvline(median_score, color='green', linestyle='--', label=f'Median: {median_score:.3f}')
                ax1.legend()
            
            # 2. Risk score by gene
            if 'Gene' in df.columns and 'score' in df.columns:
                gene_scores = df.groupby('Gene')['score'].agg(['mean', 'max', 'count']).sort_values('max', ascending=False).head(10)
                
                x_pos = np.arange(len(gene_scores))
                ax2.bar(x_pos, gene_scores['max'], alpha=0.7, color='lightcoral', label='Max Score')
                ax2.bar(x_pos, gene_scores['mean'], alpha=0.7, color='lightblue', label='Mean Score')
                
                ax2.set_xlabel('Gene')
                ax2.set_ylabel('Risk Score')
                ax2.set_title('Top 10 Genes by Risk Score')
                ax2.set_xticks(x_pos)
                ax2.set_xticklabels(gene_scores.index, rotation=45, ha='right')
                ax2.legend()
                ax2.grid(True, alpha=0.3)
            
            # 3. Risk score categories
            if 'score' in df.columns:
                # Define risk categories
                df['risk_category'] = pd.cut(df['score'], 
                                           bins=[0, 0.1, 0.3, 0.5, 0.7, 1.0], 
                                           labels=['Very Low', 'Low', 'Medium', 'High', 'Very High'])
                
                category_counts = df['risk_category'].value_counts()
                colors = ['lightgreen', 'yellow', 'orange', 'red', 'darkred']
                
                ax3.pie(category_counts.values, labels=category_counts.index, autopct='%1.1f%%', 
                       startangle=90, colors=colors[:len(category_counts)])
                ax3.set_title('Risk Score Categories')
            
            # 4. Cumulative risk distribution
            if 'score' in df.columns:
                sorted_scores = np.sort(df['score'])
                cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
                
                ax4.plot(sorted_scores, cumulative, linewidth=2, color='purple')
                ax4.set_xlabel('Risk Score')
                ax4.set_ylabel('Cumulative Probability')
                ax4.set_title('Cumulative Risk Score Distribution')
                ax4.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            output_path = self.viz_dir / f"step3_risk_distribution{vcf_basename}.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"   Risk distribution saved: {output_path}")
            return str(output_path)
            
        except Exception as e:
            print(f"   Error creating risk distribution: {e}")
            return None
    
    def _create_top_genes_analysis(self, df, vcf_basename):
        """Create top genes analysis plot"""
        try:
            if 'Gene' in df.columns and 'score' in df.columns:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
                
                # Top genes by maximum risk score
                top_genes = df.groupby('Gene')['score'].max().sort_values(ascending=False).head(15)
                
                ax1.barh(range(len(top_genes)), top_genes.values, color='lightcoral', alpha=0.7)
                ax1.set_yticks(range(len(top_genes)))
                ax1.set_yticklabels(top_genes.index)
                ax1.set_xlabel('Maximum Risk Score')
                ax1.set_title('Top 15 Genes by Risk Score')
                ax1.grid(True, alpha=0.3)
                
                # Gene association counts
                gene_counts = df['Gene'].value_counts().head(15)
                
                ax2.bar(range(len(gene_counts)), gene_counts.values, color='lightblue', alpha=0.7)
                ax2.set_xlabel('Gene')
                ax2.set_ylabel('Number of Disease Associations')
                ax2.set_title('Top 15 Genes by Association Count')
                ax2.set_xticks(range(len(gene_counts)))
                ax2.set_xticklabels(gene_counts.index, rotation=45, ha='right')
                ax2.grid(True, alpha=0.3)
                
                plt.tight_layout()
                
                output_path = self.viz_dir / f"step3_top_genes{vcf_basename}.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   Top genes analysis saved: {output_path}")
                return str(output_path)
                
        except Exception as e:
            print(f"   Error creating top genes analysis: {e}")
            return None
    
    def _create_disease_analysis(self, df, vcf_basename):
        """Create disease analysis plots"""
        try:
            if 'disease_name' in df.columns:
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
                fig.suptitle('Disease Association Analysis', fontsize=16, fontweight='bold')
                
                # 1. Top diseases by association count
                disease_counts = df['disease_name'].value_counts().head(15)
                
                ax1.barh(range(len(disease_counts)), disease_counts.values, color='lightgreen', alpha=0.7)
                ax1.set_yticks(range(len(disease_counts)))
                ax1.set_yticklabels([name[:30] + '...' if len(name) > 30 else name for name in disease_counts.index])
                ax1.set_xlabel('Number of Gene Associations')
                ax1.set_title('Top 15 Diseases by Gene Associations')
                ax1.grid(True, alpha=0.3)
                
                # 2. Disease risk score distribution
                if 'score' in df.columns:
                    disease_risk = df.groupby('disease_name')['score'].agg(['mean', 'max']).sort_values('max', ascending=False).head(10)
                    
                    x_pos = np.arange(len(disease_risk))
                    width = 0.35
                    
                    ax2.bar(x_pos - width/2, disease_risk['max'], width, label='Max Score', alpha=0.7, color='red')
                    ax2.bar(x_pos + width/2, disease_risk['mean'], width, label='Mean Score', alpha=0.7, color='blue')
                    
                    ax2.set_xlabel('Disease')
                    ax2.set_ylabel('Risk Score')
                    ax2.set_title('Top 10 Diseases by Risk Score')
                    ax2.set_xticks(x_pos)
                    ax2.set_xticklabels([name[:20] + '...' if len(name) > 20 else name for name in disease_risk.index], 
                                       rotation=45, ha='right')
                    ax2.legend()
                    ax2.grid(True, alpha=0.3)
                
                # 3. Disease category analysis (if possible to categorize)
                # Simple categorization based on keywords
                disease_categories = self._categorize_diseases(df['disease_name'].unique())
                category_counts = pd.Series(disease_categories).value_counts()
                
                ax3.pie(category_counts.values, labels=category_counts.index, autopct='%1.1f%%', startangle=90)
                ax3.set_title('Disease Categories')
                
                # 4. Genes per disease distribution
                genes_per_disease = df.groupby('disease_name')['Gene'].nunique()
                
                ax4.hist(genes_per_disease, bins=20, alpha=0.7, color='purple', edgecolor='black')
                ax4.set_xlabel('Number of Associated Genes')
                ax4.set_ylabel('Number of Diseases')
                ax4.set_title('Distribution of Genes per Disease')
                ax4.grid(True, alpha=0.3)
                
                plt.tight_layout()
                
                output_path = self.viz_dir / f"step3_disease_analysis{vcf_basename}.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   Disease analysis saved: {output_path}")
                return str(output_path)
                
        except Exception as e:
            print(f"   Error creating disease analysis: {e}")
            return None
    
    def _categorize_diseases(self, diseases):
        """Simple disease categorization based on keywords"""
        categories = []
        
        for disease in diseases:
            disease_lower = str(disease).lower()
            
            if any(keyword in disease_lower for keyword in ['cancer', 'carcinoma', 'tumor', 'malignant']):
                categories.append('Cancer')
            elif any(keyword in disease_lower for keyword in ['diabetes', 'metabolic', 'obesity']):
                categories.append('Metabolic')
            elif any(keyword in disease_lower for keyword in ['heart', 'cardiac', 'cardiovascular']):
                categories.append('Cardiovascular')
            elif any(keyword in disease_lower for keyword in ['neuro', 'brain', 'alzheimer', 'parkinson']):
                categories.append('Neurological')
            elif any(keyword in disease_lower for keyword in ['immune', 'autoimmune', 'allergy']):
                categories.append('Immune')
            elif any(keyword in disease_lower for keyword in ['bone', 'skeletal', 'muscle']):
                categories.append('Musculoskeletal')
            else:
                categories.append('Other')
        
        return categories
    
    def _create_gene_disease_heatmap(self, df, vcf_basename):
        """Create gene-disease association heatmap"""
        try:
            if 'Gene' in df.columns and 'disease_name' in df.columns and 'score' in df.columns:
                # Select top genes and diseases for heatmap
                top_genes = df.groupby('Gene')['score'].max().sort_values(ascending=False).head(15).index
                top_diseases = df.groupby('disease_name')['score'].max().sort_values(ascending=False).head(15).index
                
                # Create pivot table
                heatmap_data = df[df['Gene'].isin(top_genes) & df['disease_name'].isin(top_diseases)]
                pivot_table = heatmap_data.pivot_table(values='score', index='Gene', columns='disease_name', fill_value=0)
                
                # Truncate long disease names
                pivot_table.columns = [name[:25] + '...' if len(name) > 25 else name for name in pivot_table.columns]
                
                fig, ax = plt.subplots(1, 1, figsize=(16, 10))
                
                sns.heatmap(pivot_table, annot=True, fmt='.3f', cmap='YlOrRd', 
                           cbar_kws={'label': 'Risk Score'}, ax=ax)
                
                ax.set_title('Gene-Disease Association Heatmap (Top 15 Each)')
                ax.set_xlabel('Disease')
                ax.set_ylabel('Gene')
                
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout()
                
                output_path = self.viz_dir / f"step3_gene_disease_heatmap{vcf_basename}.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   Gene-disease heatmap saved: {output_path}")
                return str(output_path)
                
        except Exception as e:
            print(f"   Error creating gene-disease heatmap: {e}")
            return None
    
    def create_deepcv_summary(self, vcf_basename="", step_results=None):
        """
        Create comprehensive DeepCV analysis summary
        
        Args:
            vcf_basename (str): Base name for output files
            step_results (dict): Results from each step
            
        Returns:
            str: Path to summary file
        """
        print(f"\n reating DeepCV Analysis Summary")
        
        try:
            summary_file = self.viz_dir / f"deepcv_analysis_summary{vcf_basename}.txt"
            
            with open(summary_file, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("                    DeepCV ANALYSIS SUMMARY\n")
                f.write("=" * 80 + "\n\n")
                
                f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"VCF File: {vcf_basename}\n\n")
                
                # Step-by-step summary
                if step_results:
                    for step, results in step_results.items():
                        f.write(f"\n{step.upper()}:\n")
                        f.write("-" * 40 + "\n")
                        
                        if isinstance(results, dict):
                            for key, value in results.items():
                                f.write(f"  • {key}: {value}\n")
                        else:
                            f.write(f"  • {results}\n")
                
                f.write("\n" + "=" * 80 + "\n")
                f.write("                    ANALYSIS COMPLETE\n")
                f.write("=" * 80 + "\n")
            
            print(f"   Analysis summary saved: {summary_file}")
            return str(summary_file)
            
        except Exception as e:
            print(f"   Error creating analysis summary: {e}")
            return None

