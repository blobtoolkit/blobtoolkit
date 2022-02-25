
const BTK_HOST = process.env.BTK_HOST || '192.168.1.67'

const BTK_HTTPS = (String(process.env.BTK_HTTPS) === 'true')

module.exports = {
  // API URL
  'apiUrl': process.env.BTK_API_URL || (BTK_HTTPS ? 'https' : 'http') + '://'+BTK_HOST+':8000/api/v1'
}
