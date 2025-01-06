import Paper from '@mui/material/Paper';
import Text from '@mui/material/Typography';
import Link from '@mui/material/Link';
import Box from '@mui/material/Box';

const About = () => {
  return   (
  <div className="App">
  	<header className="App-header">
  		<Paper sx={{width: "100%"}}>
			<Box
				sx={{
					height: 600,
					width: '100%',
					'& .actions': {
						color: 'text.secondary',
					},
					'& .textPrimary': {
						color: 'text.primary',
					},
				}}
			>
				<h2>About</h2>
				<Text>
					<h4>References:</h4>
					[1] Glauer, Martin, et al.: Chebifier: Automating Semantic Classification in ChEBI to Accelerate
					Data-driven Discovery; Digital Discovery 3.5 (2024), <Link href="https://doi.org/10.1039/D3DD00238A">Link</Link>
					<br />
					[2] (in submission): ChemLog
					<br />
				</Text>
			</Box>
		</Paper>
    </header>
  </div>

        );
};

export default About;
