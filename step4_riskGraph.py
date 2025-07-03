#!/usr/bin/env python3
"""
DeepCV Step 4: Risk Graph Visualization with STRING API Integration
==================================================================

This module creates comprehensive visualizations including:
1. Disease-Gene-PPI networks using custom NetworkX
2. High-quality STRING network plots using STRING API
3. Comprehensive analysis reports and statistics

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import requests
import json
import time
import os
from urllib.parse import quote
from io import BytesIO
from PIL import Image
import warnings
warnings.filterwarnings('ignore')

class RiskGraphVisualizer:
    """
    Risk graph visualizer that integrates STRING database API for comprehensive network analysis
    """
    
    def __init__(self, risk_assessment_file, output_basename=""):
        """
        Initialize the risk graph visualizer
        
        Inputs:
            risk_assessment_file (str): Path to final_risk_assessment.csv
            output_basename (str): Suffix for output files
        """
        self.risk_file = risk_assessment_file
        self.output_basename = output_basename
        self.output_dir = f"results/visualizations{output_basename}"
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # STRING API configuration
        self.string_api_url = "https://string-db.org/api"
        self.species = "9606"  # Human species ID
        
        # Load and prepare data
        self.df = pd.read_csv(risk_assessment_file)
        self.top_genes = self._get_top_genes()
        
        print(f"Loaded {len(self.df)} disease associations")
        print(f"Focusing on top {len(self.top_genes)} genes")
        print(f"Output directory: {self.output_dir}")
    
    def _get_top_genes(self, n=10):
        """
        Get top N genes by disease association score
        
        Inputs:
            n (int): Number of top genes to return
            
        Outputs:
            list: Top gene names
        """
        # Get unique genes and their max scores
        gene_scores = self.df.groupby('Gene')['score'].max().sort_values(ascending=False)
        top_genes = gene_scores.head(n).index.tolist()
        
        print(f"Top {n} genes by disease association score:")
        for i, (gene, score) in enumerate(gene_scores.head(n).items(), 1):
            print(f"   {i:2d}. {gene:<12} (score: {score:.3f})")
        
        return top_genes
    
    def get_string_network_image(self, gene_list, confidence_threshold=0.7, 
                                network_type="physical", image_format="png"):
        """
        Fetch high-quality network image from STRING API
        
        Inputs:
            gene_list (list): List of gene names
            confidence_threshold (float): Minimum interaction confidence (0.0-1.0)
            network_type (str): Type of network ('physical', 'functional', 'all')
            image_format (str): Image format ('png', 'svg')
            
        Outputs:
            PIL.Image or None: Network image from STRING
        """
        # Prepare gene identifiers - STRING expects newline-separated genes
        genes_string = "%0A".join(gene_list)  # %0A is URL-encoded newline
        
        # STRING API endpoint for network images
        url = f"{self.string_api_url}/image"
        
        # Parameters for STRING network visualization
        params = {
            'identifiers': genes_string,
            'species': self.species,
            'required_score': int(confidence_threshold * 1000),  # STRING uses 0-1000 scale
            'network_flavor': network_type,
            'caller_identity': 'DeepCV_analysis'
        }
        
        try:
            print(f"Fetching STRING network for {len(gene_list)} genes...")
            print(f"   Genes: {', '.join(gene_list)}")
            print(f"   Confidence threshold: {confidence_threshold}")
            print(f"   Network type: {network_type}")
            
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200:
                # Convert response to PIL Image
                image = Image.open(BytesIO(response.content))
                print(f"STRING network image retrieved ({image.size[0]}x{image.size[1]})")
                return image
            else:
                print(f"STRING API error: {response.status_code}")
                print(f"   URL: {response.url}")
                print(f"   Response: {response.text[:200]}...")
                return None
                
        except Exception as e:
            print(f"Error fetching STRING network: {e}")
            return None
    
    def get_string_interactions(self, gene_list, confidence_threshold=0.7):
        """
        Fetch protein-protein interactions from STRING API
        
        Inputs:
            gene_list (list): List of gene names
            confidence_threshold (float): Minimum interaction confidence
            
        Outputs:
            pandas.DataFrame: Interaction data with confidence scores
        """
        # Prepare gene identifiers - STRING expects newline-separated genes
        genes_string = "%0A".join(gene_list)  # %0A is URL-encoded newline
        
        # STRING API endpoint for interactions
        url = f"{self.string_api_url}/tsv/interaction_partners"
        
        params = {
            'identifiers': genes_string,
            'species': self.species,
            'required_score': int(confidence_threshold * 1000),
            'caller_identity': 'DeepCV_analysis'
        }
        
        try:
            print(f" Fetching STRING interactions...")
            print(f"   Genes: {', '.join(gene_list)}")
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200:
                # Parse TSV response
                lines = response.text.strip().split('\n')
                if len(lines) > 1:  # Has header + data
                    # Create DataFrame from TSV
                    data = []
                    headers = lines[0].split('\t')
                    
                    for line in lines[1:]:
                        values = line.split('\t')
                        if len(values) == len(headers):
                            data.append(values)
                    
                    if data:
                        interactions_df = pd.DataFrame(data, columns=headers)
                        # Convert score to float
                        interactions_df['score'] = pd.to_numeric(interactions_df['score'], errors='coerce')
                        interactions_df = interactions_df[interactions_df['score'] >= confidence_threshold]
                        
                        print(f"Retrieved {len(interactions_df)} high-confidence interactions")
                        return interactions_df
                    
            print(f"No interactions found or API error: {response.status_code}")
            if response.status_code != 200:
                print(f"   URL: {response.url}")
                print(f"   Response: {response.text[:200]}...")
            return pd.DataFrame()
            
        except Exception as e:
            print(f" Error fetching STRING interactions: {e}")
            return pd.DataFrame()
    
    def create_string_network_plots(self, confidence_threshold=0.7):
        """
        Create enhanced network plot using STRING PPI interaction data directly
        Focus on comprehensive network only with top 10 PPIs per gene
        
        Inputs:
            confidence_threshold (float): Minimum interaction confidence
            
        Outputs:
            dict: Dictionary of plot filenames and descriptions
        """
        plots_created = {}
        
        print(f"\nSTEP 1: Creating Enhanced Network Plot from PPI Data")
        print(f" Using STRING interaction data directly (no image API)")
        print(f" Selecting top 10 PPIs per gene to reduce network complexity")
        
        # Get STRING interactions for all plots
        interactions_df = self.get_string_interactions(self.top_genes, confidence_threshold)
        
        if not interactions_df.empty:
            print(f" Retrieved {len(interactions_df)} high-confidence interactions")
            
            # Filter to top 10 PPIs per gene
            filtered_interactions_df = self._filter_top_ppis_per_gene(interactions_df, top_n=10)
            print(f" Filtered to {len(filtered_interactions_df)} top interactions ({10} per gene)")
            
            # Save PPI interactions as source-target CSV
            ppi_csv_path = os.path.join(self.output_dir, 'ppi_source_target_interactions.csv')
            self._save_ppi_source_target_csv(filtered_interactions_df, ppi_csv_path)
            
            # Create enhanced comprehensive network plot
            comprehensive_plot_path = self._create_enhanced_comprehensive_network_plot(
                filtered_interactions_df, 
                confidence_threshold
            )
            
            if comprehensive_plot_path:
                plots_created['comprehensive_network.png'] = "Enhanced Comprehensive Gene-Disease-PPI Network with Risk Scores"
                print(f"Enhanced comprehensive network saved: {comprehensive_plot_path}")
            
        else:
            print(f" No STRING interaction data available")
            print(f"   This may indicate that the genes don't have known protein interactions")
        
        return plots_created
    
    def _create_gene_info_panel(self, ax):
        """
        Create an information panel showing gene details
        
        Inputs:
            ax: Matplotlib axis for the panel
        """
        # Get gene statistics
        gene_stats = []
        for gene in self.top_genes:
            gene_data = self.df[self.df['Gene'] == gene]
            max_score = gene_data['score'].max()
            num_diseases = len(gene_data)
            top_disease = gene_data.loc[gene_data['score'].idxmax(), 'disease_name']
            
            # Truncate long disease names
            if len(top_disease) > 30:
                top_disease = top_disease[:27] + "..."
            
            gene_stats.append({
                'Gene': gene,
                'Max Score': f"{max_score:.3f}",
                'Diseases': num_diseases,
                'Top Disease': top_disease
            })
        
        # Create table
        table_data = []
        headers = ['Rank', 'Gene', 'Score', 'Diseases', 'Top Disease Association']
        
        for i, stats in enumerate(gene_stats, 1):
            table_data.append([
                f"{i:2d}",
                stats['Gene'],
                stats['Max Score'],
                f"{stats['Diseases']:2d}",
                stats['Top Disease']
            ])
        
        # Plot table
        ax.axis('tight')
        ax.axis('off')
        
        table = ax.table(cellText=table_data,
                        colLabels=headers,
                        cellLoc='left',
                        loc='center',
                        colWidths=[0.08, 0.15, 0.12, 0.12, 0.53])
        
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        
        # Style the table
        for i in range(len(headers)):
            table[(0, i)].set_facecolor('#4CAF50')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        for i in range(1, len(table_data) + 1):
            for j in range(len(headers)):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#f0f0f0')
        
        ax.set_title('Top Disease-Associated Genes', fontsize=14, fontweight='bold', pad=20)
    
    def _create_custom_string_plot(self, interactions_df, network_type, title, filename, confidence_threshold):
        """
        Create custom STRING-style network plot using interaction data
        
        Inputs:
            interactions_df (pandas.DataFrame): STRING interaction data
            network_type (str): Type of network ('physical', 'functional', 'all')
            title (str): Plot title
            filename (str): Output filename
            confidence_threshold (float): Minimum confidence threshold
            
        Outputs:
            str: Path to saved plot file
        """
        try:
            # Create network graph
            G = nx.Graph()
            
            # Add nodes for our top genes
            for gene in self.top_genes:
                G.add_node(gene, node_type='target_gene')
            
            # Add edges from interactions
            for _, row in interactions_df.iterrows():
                gene_a = row.get('preferredName_A', '')
                gene_b = row.get('preferredName_B', '')
                score = float(row.get('score', 0))
                
                if score >= confidence_threshold:
                    # Add partner genes as nodes
                    if gene_a not in G:
                        G.add_node(gene_a, node_type='partner_gene')
                    if gene_b not in G:
                        G.add_node(gene_b, node_type='partner_gene')
                    
                    # Add edge
                    G.add_edge(gene_a, gene_b, weight=score)
            
            if len(G.nodes()) == 0:
                print(f"    No network data available for {title}")
                return None
            
            # Create plot
            fig, ax = plt.subplots(1, 1, figsize=(12, 10))
            
            # Layout
            pos = nx.spring_layout(G, k=2, iterations=50)
            
            # Node colors and sizes
            node_colors = []
            node_sizes = []
            for node in G.nodes():
                if node in self.top_genes:
                    node_colors.append('#FF6B6B')  # Red for target genes
                    node_sizes.append(800)
                else:
                    node_colors.append('#4ECDC4')  # Teal for partner genes
                    node_sizes.append(400)
            
            # Draw network
            nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, 
                                 alpha=0.8, ax=ax)
            
            # Draw edges with varying thickness based on confidence
            edges = G.edges(data=True)
            edge_weights = [edge[2].get('weight', 0.5) for edge in edges]
            edge_widths = [w * 3 for w in edge_weights]  # Scale for visibility
            
            nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.6, 
                                 edge_color='gray', ax=ax)
            
            # Draw labels
            nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)
            
            # Title and formatting
            ax.set_title(f"{title}\n(Custom plot using STRING interaction data)\n"
                        f"Confidence ≥ {confidence_threshold} | {len(G.nodes())} genes | {len(G.edges())} interactions",
                        fontsize=14, fontweight='bold', pad=20)
            ax.axis('off')
            
            # Add legend
            legend_elements = [
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6B6B', 
                          markersize=10, label='Target Genes (Top Disease-Associated)'),
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#4ECDC4', 
                          markersize=8, label='Interaction Partners')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
            
            # Save plot
            output_path = os.path.join(self.output_dir, filename)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            return output_path
            
        except Exception as e:
            print(f"   Error creating custom plot: {e}")
            return None
    
    def _save_ppi_source_target_csv(self, interactions_df, output_path):
        """
        Save PPI interactions as source-target CSV file
        
        Inputs:
            interactions_df (pandas.DataFrame): STRING interaction data
            output_path (str): Path to save CSV file
        """
        try:
            # Create source-target format
            ppi_data = []
            
            for _, row in interactions_df.iterrows():
                source_gene = row.get('preferredName_A', '')
                target_gene = row.get('preferredName_B', '')
                confidence_score = float(row.get('score', 0))
                
                # Add both directions for undirected graph
                ppi_data.append({
                    'source': source_gene,
                    'target': target_gene,
                    'confidence_score': confidence_score,
                    'interaction_type': 'protein_protein_interaction'
                })
                
                # Add reverse direction
                ppi_data.append({
                    'source': target_gene,
                    'target': source_gene,
                    'confidence_score': confidence_score,
                    'interaction_type': 'protein_protein_interaction'
                })
            
            # Create DataFrame and save
            ppi_df = pd.DataFrame(ppi_data)
            ppi_df = ppi_df.drop_duplicates()  # Remove any duplicates
            ppi_df.to_csv(output_path, index=False)
            
            print(f" PPI source-target CSV saved: {output_path}")
            print(f"   Contains {len(ppi_df)} source-target pairs")
            
        except Exception as e:
            print(f" Error saving PPI CSV: {e}")
    
    def _filter_top_ppis_per_gene(self, interactions_df, top_n=10):
        """
        Filter to top N PPI interactions per target gene to reduce network complexity
        
        Inputs:
            interactions_df (pandas.DataFrame): STRING interaction data
            top_n (int): Number of top interactions to keep per gene
            
        Outputs:
            pandas.DataFrame: Filtered interaction data
        """
        try:
            filtered_interactions = []
            
            for target_gene in self.top_genes:
                # Get interactions where this gene is involved
                gene_interactions = interactions_df[
                    (interactions_df['preferredName_A'] == target_gene) | 
                    (interactions_df['preferredName_B'] == target_gene)
                ].copy()
                
                if len(gene_interactions) > 0:
                    # Sort by confidence score (descending)
                    gene_interactions = gene_interactions.sort_values('score', ascending=False)
                    
                    # Take top N interactions
                    top_interactions = gene_interactions.head(top_n)
                    filtered_interactions.append(top_interactions)
                    
                    print(f"   {target_gene}: {len(gene_interactions)} → {len(top_interactions)} interactions")
            
            if filtered_interactions:
                # Combine all filtered interactions
                result_df = pd.concat(filtered_interactions, ignore_index=True)
                # Remove duplicates (same interaction might be selected for multiple genes)
                result_df = result_df.drop_duplicates()
                return result_df
            else:
                return pd.DataFrame()
                
        except Exception as e:
            print(f" Error filtering PPIs: {e}")
            return interactions_df  # Return original if filtering fails
    
    def _create_enhanced_comprehensive_network_plot(self, interactions_df, confidence_threshold):
        """
        Create enhanced comprehensive network plot with risk scores on edges and improved layout
        
        Inputs:
            interactions_df (pandas.DataFrame): Filtered STRING interaction data
            confidence_threshold (float): Minimum confidence threshold
            
        Outputs:
            str: Path to saved enhanced comprehensive network plot
        """
        try:
            print(f"\n Creating Enhanced Comprehensive Gene-Disease-PPI Network...")
            
            # Create network graph
            G = nx.Graph()
            
            # Add disease nodes with combined risk scores
            disease_risk_scores = {}
            for _, row in self.df.iterrows():
                disease = row['disease_name']
                score = row['score']
                
                if disease in disease_risk_scores:
                    disease_risk_scores[disease] += score
                else:
                    disease_risk_scores[disease] = score
            
            # Add disease nodes
            for disease, total_score in disease_risk_scores.items():
                G.add_node(disease, 
                          node_type='disease', 
                          risk_score=total_score,
                          size=min(total_score * 1500, 1200))  # Smaller, more manageable sizes
            
            # Add gene nodes with their maximum risk scores
            gene_risk_scores = {}
            for gene in self.top_genes:
                gene_data = self.df[self.df['Gene'] == gene]
                max_score = gene_data['score'].max()
                gene_risk_scores[gene] = max_score
                
                G.add_node(gene, 
                          node_type='target_gene', 
                          risk_score=max_score,
                          size=max_score * 1000)
            
            # Add disease-gene associations with risk scores as edge weights
            for _, row in self.df.iterrows():
                gene = row['Gene']
                disease = row['disease_name']
                score = row['score']
                
                if gene in self.top_genes:
                    G.add_edge(disease, gene, 
                              edge_type='disease_association',
                              weight=score,
                              risk_score=score)
            
            # Add PPI interactions and partner genes (limited to top interactions)
            ppi_partners = set()
            for _, row in interactions_df.iterrows():
                gene_a = row.get('preferredName_A', '')
                gene_b = row.get('preferredName_B', '')
                ppi_score = float(row.get('score', 0))
                
                if ppi_score >= confidence_threshold:
                    # Add partner genes as nodes if not already present
                    if gene_a not in G:
                        G.add_node(gene_a, 
                                  node_type='ppi_partner', 
                                  risk_score=0,
                                  size=300)  # Smaller size for PPI partners
                        ppi_partners.add(gene_a)
                    
                    if gene_b not in G:
                        G.add_node(gene_b, 
                                  node_type='ppi_partner', 
                                  risk_score=0,
                                  size=300)
                        ppi_partners.add(gene_b)
                    
                    # Add PPI edge with confidence score
                    G.add_edge(gene_a, gene_b, 
                              edge_type='ppi_interaction',
                              weight=ppi_score,
                              confidence_score=ppi_score)
            
            if len(G.nodes()) == 0:
                print(f"   No network data available")
                return None
            
            # Create enhanced plot with better layout
            fig, ax = plt.subplots(1, 1, figsize=(18, 14))
            
            # Use improved hierarchical layout
            pos = self._create_improved_layout(G)
            
            # Separate nodes by type
            disease_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'disease']
            target_gene_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'target_gene']
            ppi_partner_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'ppi_partner']
            
            # Draw nodes with enhanced styling
            if disease_nodes:
                disease_sizes = [G.nodes[n].get('size', 800) for n in disease_nodes]
                nx.draw_networkx_nodes(G, pos, nodelist=disease_nodes, 
                                     node_color='#E74C3C', node_size=disease_sizes, 
                                     alpha=0.9, ax=ax, edgecolors='darkred', linewidths=1)
            
            if target_gene_nodes:
                gene_sizes = [G.nodes[n].get('size', 600) for n in target_gene_nodes]
                nx.draw_networkx_nodes(G, pos, nodelist=target_gene_nodes, 
                                     node_color='#2ECC71', node_size=gene_sizes, 
                                     alpha=0.9, ax=ax, edgecolors='darkgreen', linewidths=2)
            
            if ppi_partner_nodes:
                partner_sizes = [G.nodes[n].get('size', 300) for n in ppi_partner_nodes]
                nx.draw_networkx_nodes(G, pos, nodelist=ppi_partner_nodes, 
                                     node_color='#85C1E9', node_size=partner_sizes, 
                                     alpha=0.7, ax=ax, edgecolors='steelblue', linewidths=1)
            
            # Draw edges with risk scores and confidence scores
            disease_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'disease_association']
            ppi_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'ppi_interaction']
            
            if disease_edges:
                # Disease association edges with risk scores
                disease_weights = []
                disease_colors = []
                for u, v in disease_edges:
                    risk_score = G[u][v].get('risk_score', 0.1)
                    disease_weights.append(max(risk_score * 8, 1))  # Scale edge width by risk score
                    # Color intensity based on risk score
                    intensity = min(risk_score * 2, 1.0)
                    disease_colors.append((0.9, 0.3 - intensity * 0.2, 0.2, 0.8))  # Red with varying intensity
                
                nx.draw_networkx_edges(G, pos, edgelist=disease_edges, 
                                     width=disease_weights, edge_color=disease_colors, 
                                     ax=ax, style='-')
            
            if ppi_edges:
                # PPI edges with confidence scores
                ppi_weights = []
                ppi_colors = []
                for u, v in ppi_edges:
                    conf_score = G[u][v].get('confidence_score', 0.7)
                    ppi_weights.append(max(conf_score * 3, 0.5))  # Scale edge width by confidence
                    # Color intensity based on confidence
                    intensity = min(conf_score, 1.0)
                    ppi_colors.append((0.2, 0.7, 0.9, intensity * 0.6))  # Blue with varying intensity
                
                nx.draw_networkx_edges(G, pos, edgelist=ppi_edges, 
                                     width=ppi_weights, edge_color=ppi_colors, 
                                     ax=ax, style='-', alpha=0.6)
            
            # Draw enhanced labels
            self._draw_enhanced_labels(G, pos, ax)
            
            # Enhanced title with statistics
            ax.set_title(f"Enhanced Gene-Disease-PPI Network\n"
                        f"Risk Scores on Edges | Top 10 PPIs per Gene | Confidence ≥ {confidence_threshold}\n"
                        f"{len(disease_nodes)} Diseases | {len(target_gene_nodes)} Target Genes | {len(ppi_partner_nodes)} PPI Partners",
                        fontsize=18, fontweight='bold', pad=25)
            ax.axis('off')
            
            # Enhanced legend with risk score information
            legend_elements = [
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#E74C3C', 
                          markersize=14, label=f'Diseases ({len(disease_nodes)}) - Size ∝ Combined Risk Score'),
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ECC71', 
                          markersize=12, label=f'Target Genes ({len(target_gene_nodes)}) - Size ∝ Max Risk Score'),
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#85C1E9', 
                          markersize=8, label=f'PPI Partners ({len(ppi_partner_nodes)}) - Top 10 per Gene'),
                plt.Line2D([0], [0], color='#E74C3C', linewidth=4, 
                          label='Disease Associations - Width ∝ Risk Score'),
                plt.Line2D([0], [0], color='#85C1E9', linewidth=3, 
                          label='Protein Interactions - Width ∝ Confidence Score')
            ]
            ax.legend(handles=legend_elements, loc='upper left', fontsize=11, 
                     frameon=True, fancybox=True, shadow=True)
            
            # Enhanced statistics text box
            avg_risk = np.mean([d.get('risk_score', 0) for n, d in G.nodes(data=True) if d.get('node_type') != 'ppi_partner'])
            max_risk = max([d.get('risk_score', 0) for n, d in G.nodes(data=True)])
            
            stats_text = f"Network Statistics:\n"
            stats_text += f"• Total nodes: {len(G.nodes())}\n"
            stats_text += f"• Total edges: {len(G.edges())}\n"
            stats_text += f"• Disease associations: {len(disease_edges)}\n"
            stats_text += f"• PPI interactions: {len(ppi_edges)}\n"
            stats_text += f"• Avg. risk score: {avg_risk:.3f}\n"
            stats_text += f"• Max risk score: {max_risk:.3f}\n"
            stats_text += f"• Network density: {nx.density(G):.3f}"
            
            ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=10,
                   verticalalignment='bottom', 
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.8))
            
            # Save enhanced plot
            output_path = os.path.join(self.output_dir, 'comprehensive_network.png')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            print(f"   Enhanced network: {len(G.nodes())} nodes, {len(G.edges())} edges")
            print(f"   {len(target_gene_nodes)} target genes, {len(ppi_partner_nodes)} PPI partners")
            print(f"   {len(disease_nodes)} diseases with combined risk scores")
            print(f"   Edge widths represent risk/confidence scores")
            
            return output_path
            
        except Exception as e:
            print(f"   Error creating enhanced comprehensive network: {e}")
            return None
    
    def _create_improved_layout(self, G):
        """
        Create improved layout for better visualization
        """
        try:
            # Get nodes by type
            disease_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'disease']
            target_gene_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'target_gene']
            ppi_partner_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'ppi_partner']
            
            pos = {}
            
            # Position target genes in the center
            if target_gene_nodes:
                center_radius = 1.0
                for i, gene in enumerate(target_gene_nodes):
                    angle = 2 * np.pi * i / len(target_gene_nodes)
                    pos[gene] = (center_radius * np.cos(angle), center_radius * np.sin(angle))
            
            # Position diseases in outer ring
            if disease_nodes:
                disease_radius = 4.0
                for i, disease in enumerate(disease_nodes):
                    angle = 2 * np.pi * i / len(disease_nodes)
                    pos[disease] = (disease_radius * np.cos(angle), disease_radius * np.sin(angle))
            
            # Position PPI partners in middle ring, clustered around their target genes
            if ppi_partner_nodes:
                partner_radius = 2.5
                partners_per_gene = len(ppi_partner_nodes) // max(len(target_gene_nodes), 1)
                
                partner_idx = 0
                for gene_idx, gene in enumerate(target_gene_nodes):
                    # Find partners connected to this gene
                    gene_partners = []
                    for partner in ppi_partner_nodes:
                        if G.has_edge(gene, partner):
                            gene_partners.append(partner)
                    
                    # Position partners around their target gene
                    for j, partner in enumerate(gene_partners):
                        if partner not in pos:  # Avoid repositioning
                            base_angle = 2 * np.pi * gene_idx / len(target_gene_nodes)
                            offset_angle = (j - len(gene_partners)/2) * 0.3  # Spread around gene
                            angle = base_angle + offset_angle
                            pos[partner] = (partner_radius * np.cos(angle), partner_radius * np.sin(angle))
                
                # Position remaining partners
                for partner in ppi_partner_nodes:
                    if partner not in pos:
                        angle = 2 * np.pi * partner_idx / len(ppi_partner_nodes)
                        pos[partner] = (partner_radius * np.cos(angle), partner_radius * np.sin(angle))
                        partner_idx += 1
            
            return pos
            
        except Exception as e:
            print(f"   Layout error, using spring layout: {e}")
            return nx.spring_layout(G, k=3, iterations=50)
    
    def _draw_enhanced_labels(self, G, pos, ax):
        """
        Draw enhanced labels with better positioning and styling
        """
        try:
            # Always label target genes (most important)
            target_genes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'target_gene']
            target_labels = {gene: gene for gene in target_genes if gene in pos}
            
            # Label top 5 diseases by risk score
            disease_nodes = [(n, d.get('risk_score', 0)) for n, d in G.nodes(data=True) if d.get('node_type') == 'disease']
            disease_nodes.sort(key=lambda x: x[1], reverse=True)
            top_diseases = [n for n, _ in disease_nodes[:5]]
            disease_labels = {}
            for disease in top_diseases:
                if disease in pos:
                    # Truncate long disease names
                    label = disease[:20] + "..." if len(disease) > 20 else disease
                    disease_labels[disease] = label
            
            # Draw target gene labels (bold, larger)
            if target_labels:
                nx.draw_networkx_labels(G, pos, target_labels, font_size=12, 
                                      font_weight='bold', font_color='darkgreen', ax=ax)
            
            # Draw disease labels (smaller, different color)
            if disease_labels:
                nx.draw_networkx_labels(G, pos, disease_labels, font_size=9, 
                                      font_weight='normal', font_color='darkred', ax=ax)
            
        except Exception as e:
            print(f"   Warning: Could not draw enhanced labels: {e}")
    
    def _create_comprehensive_network_plot(self, interactions_df, confidence_threshold):
        """
        Create comprehensive network plot combining genes, diseases, and PPI with risk scores
        
        Inputs:
            interactions_df (pandas.DataFrame): STRING interaction data
            confidence_threshold (float): Minimum confidence threshold
            
        Outputs:
            str: Path to saved comprehensive network plot
        """
        try:
            print(f"\n Creating Comprehensive Gene-Disease-PPI Network...")
            
            # Create network graph
            G = nx.Graph()
            
            # Add disease nodes with combined risk scores
            disease_risk_scores = {}
            for _, row in self.df.iterrows():
                disease = row['disease_name']
                score = row['score']
                
                if disease in disease_risk_scores:
                    disease_risk_scores[disease] += score
                else:
                    disease_risk_scores[disease] = score
            
            # Add disease nodes
            for disease, total_score in disease_risk_scores.items():
                G.add_node(disease, 
                          node_type='disease', 
                          risk_score=total_score,
                          size=min(total_score * 3000, 2000))  # Scale size but cap it
            
            # Add gene nodes with their maximum risk scores
            gene_risk_scores = {}
            for gene in self.top_genes:
                gene_data = self.df[self.df['Gene'] == gene]
                max_score = gene_data['score'].max()
                gene_risk_scores[gene] = max_score
                
                G.add_node(gene, 
                          node_type='target_gene', 
                          risk_score=max_score,
                          size=max_score * 2000)
            
            # Add disease-gene associations
            for _, row in self.df.iterrows():
                gene = row['Gene']
                disease = row['disease_name']
                score = row['score']
                
                if gene in self.top_genes:
                    G.add_edge(disease, gene, 
                              edge_type='disease_association',
                              weight=score)
            
            # Add PPI interactions and partner genes
            ppi_partners = set()
            for _, row in interactions_df.iterrows():
                gene_a = row.get('preferredName_A', '')
                gene_b = row.get('preferredName_B', '')
                ppi_score = float(row.get('score', 0))
                
                if ppi_score >= confidence_threshold:
                    # Add partner genes as nodes if not already present
                    if gene_a not in G:
                        G.add_node(gene_a, 
                                  node_type='ppi_partner', 
                                  risk_score=0,
                                  size=400)
                        ppi_partners.add(gene_a)
                    
                    if gene_b not in G:
                        G.add_node(gene_b, 
                                  node_type='ppi_partner', 
                                  risk_score=0,
                                  size=400)
                        ppi_partners.add(gene_b)
                    
                    # Add PPI edge
                    G.add_edge(gene_a, gene_b, 
                              edge_type='ppi_interaction',
                              weight=ppi_score)
            
            if len(G.nodes()) == 0:
                print(f"   No network data available")
                return None
            
            # Create comprehensive plot
            fig, ax = plt.subplots(1, 1, figsize=(16, 12))
            
            # Use hierarchical layout
            pos = self._create_hierarchical_layout(G)
            
            # Separate nodes by type
            disease_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'disease']
            target_gene_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'target_gene']
            ppi_partner_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'ppi_partner']
            
            # Draw nodes with different colors and sizes
            if disease_nodes:
                disease_sizes = [G.nodes[n].get('size', 1000) for n in disease_nodes]
                nx.draw_networkx_nodes(G, pos, nodelist=disease_nodes, 
                                     node_color='#FF6B6B', node_size=disease_sizes, 
                                     alpha=0.8, ax=ax, label='Diseases')
            
            if target_gene_nodes:
                gene_sizes = [G.nodes[n].get('size', 800) for n in target_gene_nodes]
                nx.draw_networkx_nodes(G, pos, nodelist=target_gene_nodes, 
                                     node_color='#4ECDC4', node_size=gene_sizes, 
                                     alpha=0.8, ax=ax, label='Target Genes')
            
            if ppi_partner_nodes:
                partner_sizes = [G.nodes[n].get('size', 400) for n in ppi_partner_nodes]
                nx.draw_networkx_nodes(G, pos, nodelist=ppi_partner_nodes, 
                                     node_color='#95E1D3', node_size=partner_sizes, 
                                     alpha=0.6, ax=ax, label='PPI Partners')
            
            # Draw edges with different colors
            disease_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'disease_association']
            ppi_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'ppi_interaction']
            
            if disease_edges:
                disease_weights = [G[u][v].get('weight', 0.5) * 5 for u, v in disease_edges]
                nx.draw_networkx_edges(G, pos, edgelist=disease_edges, 
                                     width=disease_weights, edge_color='#FF6B6B', 
                                     alpha=0.7, ax=ax)
            
            if ppi_edges:
                ppi_weights = [G[u][v].get('weight', 0.5) * 3 for u, v in ppi_edges]
                nx.draw_networkx_edges(G, pos, edgelist=ppi_edges, 
                                     width=ppi_weights, edge_color='#4ECDC4', 
                                     alpha=0.5, ax=ax)
            
            # Draw labels with smart positioning
            self._draw_smart_labels(G, pos, ax)
            
            # Title and formatting
            ax.set_title(f"Comprehensive Gene-Disease-PPI Network\n"
                        f"Combined Risk Scores | Confidence ≥ {confidence_threshold}\n"
                        f"{len(disease_nodes)} Diseases | {len(target_gene_nodes)} Target Genes | {len(ppi_partner_nodes)} PPI Partners",
                        fontsize=16, fontweight='bold', pad=20)
            ax.axis('off')
            
            # Add comprehensive legend
            legend_elements = [
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6B6B', 
                          markersize=12, label=f'Diseases ({len(disease_nodes)})'),
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#4ECDC4', 
                          markersize=10, label=f'Target Genes ({len(target_gene_nodes)})'),
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#95E1D3', 
                          markersize=8, label=f'PPI Partners ({len(ppi_partner_nodes)})'),
                plt.Line2D([0], [0], color='#FF6B6B', linewidth=3, 
                          label='Disease Associations'),
                plt.Line2D([0], [0], color='#4ECDC4', linewidth=2, 
                          label='Protein Interactions')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
            
            # Add statistics text box
            stats_text = f"Network Statistics:\n"
            stats_text += f"• Total nodes: {len(G.nodes())}\n"
            stats_text += f"• Total edges: {len(G.edges())}\n"
            stats_text += f"• Disease associations: {len(disease_edges)}\n"
            stats_text += f"• PPI interactions: {len(ppi_edges)}\n"
            stats_text += f"• Avg. risk score: {np.mean([d.get('risk_score', 0) for n, d in G.nodes(data=True)]):.3f}"
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            # Save plot
            output_path = os.path.join(self.output_dir, 'comprehensive_network.png')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"   Network contains {len(G.nodes())} nodes and {len(G.edges())} edges")
            print(f"   {len(target_gene_nodes)} target genes, {len(ppi_partner_nodes)} PPI partners")
            print(f"   {len(disease_nodes)} diseases with combined risk scores")
            
            return output_path
            
        except Exception as e:
            print(f"   Error creating comprehensive network: {e}")
            return None
    
    def _create_hierarchical_layout(self, G):
        """
        Create hierarchical layout for the comprehensive network
        """
        try:
            # Try to use hierarchical layout if possible
            pos = {}
            
            # Get nodes by type
            disease_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'disease']
            target_gene_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'target_gene']
            ppi_partner_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'ppi_partner']
            
            # Position diseases at the top
            for i, disease in enumerate(disease_nodes):
                angle = 2 * np.pi * i / len(disease_nodes) if len(disease_nodes) > 1 else 0
                pos[disease] = (3 * np.cos(angle), 3 + 2 * np.sin(angle))
            
            # Position target genes in the middle
            for i, gene in enumerate(target_gene_nodes):
                angle = 2 * np.pi * i / len(target_gene_nodes) if len(target_gene_nodes) > 1 else 0
                pos[gene] = (1.5 * np.cos(angle), 1.5 * np.sin(angle))
            
            # Position PPI partners at the bottom
            for i, partner in enumerate(ppi_partner_nodes):
                angle = 2 * np.pi * i / len(ppi_partner_nodes) if len(ppi_partner_nodes) > 1 else 0
                pos[partner] = (2 * np.cos(angle), -2 + 1.5 * np.sin(angle))
            
            return pos
            
        except:
            # Fallback to spring layout
            return nx.spring_layout(G, k=3, iterations=50)
    
    def _draw_smart_labels(self, G, pos, ax):
        """
        Draw labels with smart positioning to avoid overlap
        """
        try:
            # Draw labels for important nodes only
            important_nodes = []
            
            # Always label target genes
            target_genes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'target_gene']
            important_nodes.extend(target_genes)
            
            # Label top diseases by risk score
            disease_nodes = [(n, d.get('risk_score', 0)) for n, d in G.nodes(data=True) if d.get('node_type') == 'disease']
            disease_nodes.sort(key=lambda x: x[1], reverse=True)
            top_diseases = [n for n, _ in disease_nodes[:5]]  # Top 5 diseases
            important_nodes.extend(top_diseases)
            
            # Create label dict for important nodes only
            labels = {node: node for node in important_nodes if node in pos}
            
            # Truncate long labels
            for node in labels:
                if len(labels[node]) > 15:
                    labels[node] = labels[node][:12] + "..."
            
            nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold', ax=ax)
            
        except Exception as e:
            print(f"   Warning: Could not draw labels: {e}")
    
    def create_disease_gene_network(self, confidence_threshold=0.7):
        """
        Create enhanced Disease-Gene-PPI network visualization with improved layout and labels
        
        Inputs:
            confidence_threshold (float): Minimum PPI confidence
            
        Outputs:
            str: Path to saved network plot
        """
        print(f"\n Creating Enhanced Disease-Gene-PPI Network...")
        
        # Get STRING interactions for top genes
        interactions_df = self.get_string_interactions(self.top_genes, confidence_threshold)
        
        # Create network graph
        G = nx.Graph()
        
        # Add disease nodes with combined risk scores
        disease_risk_scores = {}
        for _, row in self.df.iterrows():
            disease = row['disease_name']
            score = row['score']
            if disease in disease_risk_scores:
                disease_risk_scores[disease] += score
            else:
                disease_risk_scores[disease] = score
        
        for disease, total_score in disease_risk_scores.items():
            G.add_node(disease, node_type='disease', 
                      size=min(total_score * 1200, 1000),  # Cap the size
                      risk_score=total_score)
        
        # Add gene nodes
        for gene in self.top_genes:
            gene_data = self.df[self.df['Gene'] == gene]
            max_score = gene_data['score'].max()
            G.add_node(gene, node_type='gene', 
                      size=max_score * 1500,
                      risk_score=max_score)
        
        # Add disease-gene edges
        for _, row in self.df.iterrows():
            if row['Gene'] in self.top_genes:
                G.add_edge(row['disease_name'], row['Gene'], 
                          edge_type='disease_association', 
                          weight=row['score'])
        
        # Add PPI edges from STRING (only between target genes)
        if not interactions_df.empty:
            for _, interaction in interactions_df.iterrows():
                gene1 = interaction.get('preferredName_A', interaction.get('stringId_A', ''))
                gene2 = interaction.get('preferredName_B', interaction.get('stringId_B', ''))
                score = float(interaction.get('score', 0))
                
                if gene1 in self.top_genes and gene2 in self.top_genes:
                    G.add_edge(gene1, gene2, 
                              edge_type='ppi', 
                              weight=score)
        
        # Create enhanced visualization
        fig, ax = plt.subplots(1, 1, figsize=(18, 14))
        
        # Use improved layout with better separation
        pos = self._create_enhanced_disease_gene_layout(G)
        
        # Separate nodes by type
        disease_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'disease']
        gene_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'gene']
        
        # Draw disease nodes with enhanced styling
        if disease_nodes:
            disease_sizes = [G.nodes[n].get('size', 800) for n in disease_nodes]
            nx.draw_networkx_nodes(G, pos, nodelist=disease_nodes, 
                                  node_color='#E74C3C', node_size=disease_sizes, 
                                  alpha=0.8, ax=ax, edgecolors='darkred', linewidths=1.5)
        
        # Draw gene nodes with enhanced styling
        if gene_nodes:
            gene_sizes = [G.nodes[n].get('size', 600) for n in gene_nodes]
            nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, 
                                  node_color='#2ECC71', node_size=gene_sizes, 
                                  alpha=0.9, ax=ax, edgecolors='darkgreen', linewidths=2)
        
        # Draw edges with enhanced styling
        disease_edges = [(u, v) for u, v, d in G.edges(data=True) 
                        if d.get('edge_type') == 'disease_association']
        ppi_edges = [(u, v) for u, v, d in G.edges(data=True) 
                    if d.get('edge_type') == 'ppi']
        
        if disease_edges:
            # Disease association edges with varying widths based on risk scores
            disease_weights = [G[u][v].get('weight', 0.1) * 6 for u, v in disease_edges]
            nx.draw_networkx_edges(G, pos, edgelist=disease_edges, 
                                  edge_color='#E74C3C', alpha=0.7, 
                                  width=disease_weights, ax=ax)
        
        if ppi_edges:
            # PPI edges with varying widths based on confidence
            ppi_weights = [G[u][v].get('weight', 0.7) * 4 for u, v in ppi_edges]
            nx.draw_networkx_edges(G, pos, edgelist=ppi_edges, 
                                  edge_color='#2ECC71', alpha=0.8, 
                                  width=ppi_weights, ax=ax)
        
        # Draw enhanced labels with better positioning
        self._draw_enhanced_disease_gene_labels(G, pos, ax)
        
        # Enhanced title
        ax.set_title(f"Enhanced Disease-Gene-PPI Network\n"
                    f"Top {len(self.top_genes)} Genes with High-Confidence Interactions (≥{confidence_threshold})\n"
                    f"{len(disease_nodes)} Diseases | {len(gene_nodes)} Genes | {len(ppi_edges)} PPI Interactions",
                    fontsize=16, fontweight='bold', pad=25)
        
        # Enhanced legend positioned to avoid overlap
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#E74C3C', 
                      markersize=12, label=f'Diseases ({len(disease_nodes)})'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ECC71', 
                      markersize=10, label=f'Target Genes ({len(gene_nodes)})'),
            plt.Line2D([0], [0], color='#E74C3C', linewidth=3, 
                      label='Disease Associations'),
            plt.Line2D([0], [0], color='#2ECC71', linewidth=3, 
                      label='Protein Interactions')
        ]
        
        # Position legend in bottom left to avoid overlap
        ax.legend(handles=legend_elements, loc='lower left', fontsize=11, 
                 frameon=True, fancybox=True, shadow=True, 
                 bbox_to_anchor=(0.02, 0.02))
        
        # Add statistics box in top left
        stats_text = f"Network Statistics:\n"
        stats_text += f"• Total nodes: {len(G.nodes())}\n"
        stats_text += f"• Disease associations: {len(disease_edges)}\n"
        stats_text += f"• PPI interactions: {len(ppi_edges)}\n"
        stats_text += f"• Network density: {nx.density(G):.3f}"
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
               verticalalignment='top', 
               bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.8))
        
        ax.axis('off')
        
        # Save enhanced plot
        output_path = os.path.join(self.output_dir, 'disease_gene_ppi_network.png')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f" Enhanced Disease-Gene-PPI network saved: {output_path}")
        print(f"   {len(disease_nodes)} diseases, {len(gene_nodes)} genes")
        print(f"   {len(disease_edges)} disease associations, {len(ppi_edges)} PPI interactions")
        
        return output_path
    
    def _create_enhanced_disease_gene_layout(self, G):
        """
        Create enhanced layout for disease-gene network with better separation
        """
        try:
            # Get nodes by type
            disease_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'disease']
            gene_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'gene']
            
            pos = {}
            
            # Position genes in the center in a circle
            if gene_nodes:
                center_radius = 1.5
                for i, gene in enumerate(gene_nodes):
                    angle = 2 * np.pi * i / len(gene_nodes)
                    pos[gene] = (center_radius * np.cos(angle), center_radius * np.sin(angle))
            
            # Position diseases in outer rings, grouped by connection to genes
            if disease_nodes:
                disease_radius = 4.5
                
                # Group diseases by their connected genes
                disease_gene_connections = {}
                for disease in disease_nodes:
                    connected_genes = []
                    for gene in gene_nodes:
                        if G.has_edge(disease, gene):
                            connected_genes.append(gene)
                    disease_gene_connections[disease] = connected_genes
                
                # Position diseases around their connected genes
                positioned_diseases = set()
                for gene_idx, gene in enumerate(gene_nodes):
                    # Find diseases connected to this gene
                    connected_diseases = [d for d, genes in disease_gene_connections.items() 
                                        if gene in genes and d not in positioned_diseases]
                    
                    if connected_diseases:
                        base_angle = 2 * np.pi * gene_idx / len(gene_nodes)
                        
                        for j, disease in enumerate(connected_diseases):
                            # Spread diseases around the gene
                            offset_angle = (j - len(connected_diseases)/2) * 0.4
                            angle = base_angle + offset_angle
                            pos[disease] = (disease_radius * np.cos(angle), 
                                          disease_radius * np.sin(angle))
                            positioned_diseases.add(disease)
                
                # Position remaining diseases
                remaining_diseases = [d for d in disease_nodes if d not in positioned_diseases]
                for i, disease in enumerate(remaining_diseases):
                    angle = 2 * np.pi * i / len(remaining_diseases) + np.pi/4  # Offset to avoid overlap
                    pos[disease] = (disease_radius * np.cos(angle), disease_radius * np.sin(angle))
            
            return pos
            
        except Exception as e:
            print(f"   Layout error, using spring layout: {e}")
            return nx.spring_layout(G, k=4, iterations=100)
    
    def _draw_enhanced_disease_gene_labels(self, G, pos, ax):
        """
        Draw enhanced labels with smart positioning to avoid overlap
        """
        try:
            # Always label genes (most important)
            gene_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'gene']
            gene_labels = {gene: gene for gene in gene_nodes if gene in pos}
            
            # Label only top diseases by risk score to avoid clutter
            disease_nodes = [(n, d.get('risk_score', 0)) for n, d in G.nodes(data=True) 
                          if d.get('node_type') == 'disease']
            disease_nodes.sort(key=lambda x: x[1], reverse=True)
            
            # Select top diseases but limit to avoid overcrowding
            max_disease_labels = min(8, len(disease_nodes))
            top_diseases = [n for n, _ in disease_nodes[:max_disease_labels]]
            
            disease_labels = {}
            for disease in top_diseases:
                if disease in pos:
                    # Truncate long disease names more aggressively
                    if len(disease) > 25:
                        label = disease[:22] + "..."
                    else:
                        label = disease
                    disease_labels[disease] = label
            
            # Draw gene labels (bold, larger, green)
            if gene_labels:
                nx.draw_networkx_labels(G, pos, gene_labels, font_size=11, 
                                      font_weight='bold', font_color='darkgreen', ax=ax)
            
            # Draw disease labels (smaller, red) with offset positioning to reduce overlap
            if disease_labels:
                # Create offset positions for disease labels
                offset_pos = {}
                for disease in disease_labels:
                    x, y = pos[disease]
                    # Offset labels slightly outward from center
                    offset_factor = 1.15
                    offset_pos[disease] = (x * offset_factor, y * offset_factor)
                
                nx.draw_networkx_labels(G, offset_pos, disease_labels, font_size=8, 
                                      font_weight='normal', font_color='darkred', ax=ax)
            
        except Exception as e:
            print(f"   ⚠️ Warning: Could not draw enhanced labels: {e}")

    
    def create_analysis_summary(self):
        """
        Create comprehensive analysis summary with statistics
        
        Outputs:
            dict: Summary statistics and file paths
        """
        print(f"\n📊 Creating analysis summary...")
        
        # Calculate statistics
        stats = {
            'total_associations': len(self.df),
            'unique_genes': self.df['Gene'].nunique(),
            'unique_diseases': self.df['disease_name'].nunique(),
            'avg_score': self.df['score'].mean(),
            'max_score': self.df['score'].max(),
            'top_gene': self.df.loc[self.df['score'].idxmax(), 'Gene'],
            'top_disease': self.df.loc[self.df['score'].idxmax(), 'disease_name']
        }
        
        # Create summary plot
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # Top diseases by gene count
        disease_counts = self.df['disease_name'].value_counts().head(10)
        ax1.barh(range(len(disease_counts)), disease_counts.values, color='lightcoral')
        ax1.set_yticks(range(len(disease_counts)))
        ax1.set_yticklabels([d[:30] + "..." if len(d) > 30 else d for d in disease_counts.index])
        ax1.set_xlabel('Number of Associated Genes')
        ax1.set_title('Top Diseases by Gene Count')
        
        # Score distribution
        ax2.hist(self.df['score'], bins=20, color='skyblue', alpha=0.7, edgecolor='black')
        ax2.set_xlabel('Disease Association Score')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Distribution of Association Scores')
        ax2.axvline(stats['avg_score'], color='red', linestyle='--', label=f'Mean: {stats["avg_score"]:.3f}')
        ax2.legend()
        
        # Top genes by max score
        gene_max_scores = self.df.groupby('Gene')['score'].max().sort_values(ascending=False).head(10)
        ax3.bar(range(len(gene_max_scores)), gene_max_scores.values, color='lightgreen')
        ax3.set_xticks(range(len(gene_max_scores)))
        ax3.set_xticklabels(gene_max_scores.index, rotation=45, ha='right')
        ax3.set_ylabel('Maximum Association Score')
        ax3.set_title('Top Genes by Association Score')
        
        # Statistics text
        ax4.axis('off')
        stats_text = f"""
        📊 ANALYSIS SUMMARY
        ==================
        
        Total Disease Associations: {stats['total_associations']:,}
        Unique Genes: {stats['unique_genes']:,}
        Unique Diseases: {stats['unique_diseases']:,}
        
        Average Association Score: {stats['avg_score']:.3f}
        Maximum Association Score: {stats['max_score']:.3f}
        
        Top Gene: {stats['top_gene']}
        Top Disease: {stats['top_disease'][:40]}...
        
        Analysis completed successfully!
        High-confidence interactions retrieved from STRING database.
        """
        ax4.text(0.1, 0.9, stats_text, transform=ax4.transAxes, fontsize=12,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.suptitle('DeepCV Analysis Summary', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # Save summary plot
        summary_path = os.path.join(self.output_dir, 'analysis_summary.png')
        plt.savefig(summary_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✅ Analysis summary saved: {summary_path}")
        return stats
    
    def save_detailed_reports(self, interactions_df):
        """
        Save detailed CSV and text reports
        
        Inputs:
            interactions_df (pandas.DataFrame): STRING interaction data
        """
        print(f"\n📄 Saving detailed reports...")
        
        # Gene summary
        gene_summary = self.df.groupby('Gene').agg({
            'score': ['count', 'mean', 'max'],
            'disease_name': 'nunique'
        }).round(3)
        gene_summary.columns = ['Num_Associations', 'Avg_Score', 'Max_Score', 'Unique_Diseases']
        gene_summary = gene_summary.sort_values('Max_Score', ascending=False)
        gene_summary_path = os.path.join(self.output_dir, 'gene_summary.csv')
        gene_summary.to_csv(gene_summary_path)
        
        # Disease summary
        disease_summary = self.df.groupby('disease_name').agg({
            'score': ['count', 'mean', 'max'],
            'Gene': 'nunique'
        }).round(3)
        disease_summary.columns = ['Num_Associations', 'Avg_Score', 'Max_Score', 'Unique_Genes']
        disease_summary = disease_summary.sort_values('Max_Score', ascending=False)
        disease_summary_path = os.path.join(self.output_dir, 'disease_summary.csv')
        disease_summary.to_csv(disease_summary_path)
        
        # High-confidence associations
        high_conf = self.df[self.df['score'] >= 0.3].copy()
        high_conf_path = os.path.join(self.output_dir, 'high_confidence_associations.csv')
        high_conf.to_csv(high_conf_path, index=False)
        
        # STRING interactions
        if not interactions_df.empty:
            interactions_path = os.path.join(self.output_dir, 'string_interactions.csv')
            interactions_df.to_csv(interactions_path, index=False)
            print(f"✅ STRING interactions saved: {interactions_path}")
        
        print(f"✅ Gene summary saved: {gene_summary_path}")
        print(f"✅ Disease summary saved: {disease_summary_path}")
        print(f"✅ High-confidence associations saved: {high_conf_path}")
    
    def run_complete_analysis(self, confidence_threshold=0.7):
        """
        Run complete risk graph analysis with STRING API integration
        
        Inputs:
            confidence_threshold (float): Minimum interaction confidence
            
        Outputs:
            dict: Analysis results and file paths
        """
        print(f"\n🚀 Starting Complete DeepCV Risk Graph Analysis")
        print(f"=" * 70)
        
        results = {
            'string_plots': {},
            'custom_plots': [],
            'reports': [],
            'statistics': {}
        }
        
        # 1. Create STRING network plots
        print(f"\n🌐 STEP 1: Creating STRING Network Plots")
        results['string_plots'] = self.create_string_network_plots(confidence_threshold)
        
        # 2. Get STRING interactions for custom analysis
        print(f"\n🔗 STEP 2: Fetching STRING Interactions")
        interactions_df = self.get_string_interactions(self.top_genes, confidence_threshold)
        
        # 3. Create custom Disease-Gene-PPI network
        print(f"\n🎨 STEP 3: Creating Custom Disease-Gene-PPI Network")
        custom_network_path = self.create_disease_gene_network(confidence_threshold)
        results['custom_plots'].append(custom_network_path)
        
        # 4. Create analysis summary
        print(f"\n📊 STEP 4: Creating Analysis Summary")
        results['statistics'] = self.create_analysis_summary()
        
        # 5. Save detailed reports
        print(f"\n📄 STEP 5: Saving Detailed Reports")
        self.save_detailed_reports(interactions_df)
        
        # Final summary
        print(f"\n" + "=" * 70)
        print(f"🎉 RISK GRAPH ANALYSIS COMPLETE!")
        print(f"=" * 70)
        print(f"📁 Output directory: {self.output_dir}")
        print(f"🌐 STRING plots created: {len(results['string_plots'])}")
        print(f"🎨 Custom plots created: {len(results['custom_plots'])}")
        print(f"📊 Total associations analyzed: {results['statistics']['total_associations']}")
        print(f"🧬 Unique genes: {results['statistics']['unique_genes']}")
        print(f"🏥 Unique diseases: {results['statistics']['unique_diseases']}")
        print(f"=" * 70)
        
        return results

def main(risk_assessment_file=None, output_basename="", confidence_threshold=0.7):
    """
    Main function to run risk graph visualization analysis
    
    Inputs:
        risk_assessment_file (str): Path to final_risk_assessment.csv
        output_basename (str): Suffix for output files
        confidence_threshold (float): Minimum interaction confidence
    """
    # Default file path if not provided
    if risk_assessment_file is None:
        risk_assessment_file = "results/final_risk_assessment.csv"
        # Look for files with basename suffix
        import glob
        pattern = f"results/final_risk_assessment*{output_basename}.csv"
        matching_files = glob.glob(pattern)
        if matching_files:
            risk_assessment_file = matching_files[0]
    
    # Check if file exists
    if not os.path.exists(risk_assessment_file):
        print(f"❌ Error: Risk assessment file not found: {risk_assessment_file}")
        print(f"   Please run DeepCV steps 1-3 first to generate the risk assessment file.")
        return None
    
    # Initialize visualizer
    visualizer = RiskGraphVisualizer(risk_assessment_file, output_basename)
    
    # Run complete analysis
    results = visualizer.run_complete_analysis(confidence_threshold)
    
    return results

if __name__ == "__main__":
    import sys
    
    # Command line usage
    if len(sys.argv) > 1:
        risk_file = sys.argv[1]
        basename = sys.argv[2] if len(sys.argv) > 2 else ""
        confidence = float(sys.argv[3]) if len(sys.argv) > 3 else 0.7
        main(risk_file, basename, confidence)
    else:
        # Auto-detect files
        main()

